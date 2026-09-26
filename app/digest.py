"""Core restriction-digest logic (no web framework dependencies)."""

from __future__ import annotations

from dataclasses import dataclass

# IUPAC degenerate base codes -> the concrete bases they match.
IUPAC_MATCH: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("GC"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}

# Complement of each IUPAC code (used to build the reverse complement site).
IUPAC_COMPLEMENT: dict[str, str] = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "V": "B",
    "D": "H",
    "H": "D",
    "N": "N",
}


@dataclass(frozen=True)
class Fragment:
    """A digest product: zero-based start on the original sequence, length,
    and the bases in the original (5'->3') direction."""

    start: int
    length: int
    sequence: str


def reverse_complement(site: str) -> str:
    """Reverse complement of an IUPAC recognition string."""
    return "".join(IUPAC_COMPLEMENT[base] for base in reversed(site))


def _match_starts(sequence: str, pattern: str, circular: bool) -> list[int]:
    """All start positions where ``pattern`` matches ``sequence``.

    Overlapping matches are kept.  For circular molecules a match may span
    the origin, so starts run over every position and the tail of the
    sequence is prepended conceptually by scanning an extended string.
    """
    n = len(sequence)
    m = len(pattern)
    if circular:
        haystack = sequence + sequence[: m - 1]
        limit = n
    else:
        haystack = sequence
        limit = n - m + 1

    starts: list[int] = []
    for s in range(limit):
        for i, code in enumerate(pattern):
            if haystack[s + i] not in IUPAC_MATCH[code]:
                break
        else:
            starts.append(s)
    return starts


def enzyme_cut_positions(
    sequence: str, site: str, cut_offset: int, circular: bool
) -> list[int]:
    """Sorted, de-duplicated zero-based cut positions for one enzyme.

    The recognition string and its reverse complement are both searched:
    a forward-strand match at start ``s`` cuts at ``s + k``; a
    reverse-complement match cuts at ``s + len(site) - k``.  On circular
    molecules positions are taken modulo the sequence length.
    """
    n = len(sequence)
    m = len(site)
    cuts: set[int] = set()

    for s in _match_starts(sequence, site, circular):
        cuts.add((s + cut_offset) % n if circular else s + cut_offset)

    rc_site = reverse_complement(site)
    for s in _match_starts(sequence, rc_site, circular):
        cuts.add((s + m - cut_offset) % n if circular else s + m - cut_offset)

    return sorted(cuts)


def digest(sequence: str, cut_positions: list[int], circular: bool) -> list[Fragment]:
    """Cut ``sequence`` at all ``cut_positions`` simultaneously.

    Fragments are returned sorted by start position, each carrying its
    bases in the original sequence direction.  With no cut positions the
    original sequence is returned as a single fragment; empty fragments
    are never produced.
    """
    n = len(sequence)
    cuts = sorted(set(cut_positions))

    if circular:
        cuts = sorted({c % n for c in cuts})
        if not cuts:
            return [Fragment(start=0, length=n, sequence=sequence)]
        fragments: list[Fragment] = []
        for i, c in enumerate(cuts):
            nxt = cuts[(i + 1) % len(cuts)]
            length = (nxt - c) % n
            if length == 0:
                # A single cut linearises the circle at that point.
                length = n
            if nxt <= c:  # wraps around the origin
                bases = sequence[c:] + sequence[:nxt]
            else:
                bases = sequence[c:nxt]
            fragments.append(Fragment(start=c, length=length, sequence=bases))
        return fragments

    boundaries = sorted({0, n, *cuts})
    return [
        Fragment(start=a, length=b - a, sequence=sequence[a:b])
        for a, b in zip(boundaries, boundaries[1:])
        if b > a
    ]


def cut_probabilities(
    positions_by_enzyme: list[tuple[list[int], float]],
) -> dict[int, float]:
    """Per-position effective cut probability when every enzyme cuts
    independently at each of its (already de-duplicated) positions.

    Each entry ``(positions, p)`` describes one enzyme: every listed
    position is independently cut with probability ``p``.  Positions hit
    by several enzymes form a cut as soon as any of those enzymes cuts,
    so the combined probability is ``1 - prod(1 - p)``.
    """
    stay_alive: dict[int, float] = {}
    for positions, p in positions_by_enzyme:
        if p == 0.0:
            continue
        survive = 1.0 - p
        for pos in positions:
            stay_alive[pos] = stay_alive.get(pos, 1.0) * survive
    return {pos: 1.0 - stay for pos, stay in stay_alive.items() if stay < 1.0}


def _prefix_suffix_products(values: list[float]) -> tuple[list[float], list[float]]:
    """Prefix/suffix product arrays, division-free (values may contain 0)."""
    prefix = [1.0] * (len(values) + 1)
    for i, v in enumerate(values):
        prefix[i + 1] = prefix[i] * v
    suffix = [1.0] * (len(values) + 1)
    for i in range(len(values) - 1, -1, -1):
        suffix[i] = suffix[i + 1] * values[i]
    return prefix, suffix


def partial_digest_distribution(
    n: int, cuts: dict[int, float], circular: bool
) -> tuple[dict[int, float], float]:
    """Expected fragment-length distribution of a probabilistic digest.

    ``cuts`` maps each candidate cut position to its probability ``q``
    (0 < q <= 1); cut events at distinct positions are independent.
    Returns ``(length -> expected number of fragments of that length,
    probability that the molecule has no cut at all)``.

    Only positive-expectation lengths are present.  For a circle, both
    the no-cut outcome and the exactly-one-cut outcome yield one
    full-length fragment, while the no-cut probability only counts the
    former.

    Given K candidate positions the distribution itself is computed in
    O(K) memory and O(K^2) time without sampling or enumerating cut
    subsets: every expectation is a sum of products over the outcomes
    of independent Bernoulli variables, accumulated by length.
    """
    if circular:
        # Residues modulo n identify the same physical position, so
        # probabilities landing on the same residue are combined as
        # independent cut events.
        merged: dict[int, float] = {}
        for pos, q in cuts.items():
            residue = pos % n
            merged[residue] = 1.0 - (1.0 - merged.get(residue, 0.0)) * (1.0 - q)
        cuts = {p: q for p, q in merged.items() if 0.0 < q <= 1.0}

    positions = sorted(cuts)
    probabilities = [cuts[p] for p in positions]
    stays = [1.0 - q for q in probabilities]

    expected: dict[int, float] = {}

    def add(length: int, weight: float) -> None:
        if weight:
            expected[length] = expected.get(length, 0.0) + weight

    if not circular:
        # The molecular ends (0 and n) are never random cut points.
        inner = [(p, q) for p, q in zip(positions, probabilities) if 0 < p < n]
        if not inner:
            return {n: 1.0}, 1.0

        cut_positions = [p for p, _ in inner]
        cut_probs = [q for _, q in inner]
        cut_stays = [1.0 - q for q in cut_probs]
        c_prefix, c_suffix = _prefix_suffix_products(cut_stays)
        no_cut_prob = c_prefix[-1]

        # Prefix fragments [0, p_i): length p_i.
        for i, (p, q) in enumerate(inner):
            add(p, q * c_prefix[i])

        # Suffix fragments [p_i, n): length n - p_i.
        for i, (p, q) in enumerate(inner):
            add(n - p, q * c_suffix[i + 1])

        # Interior fragments [p_i, p_j), i < j.  The probability that a
        # fixed pair of positions are consecutive cuts equals the
        # product of stay-probabilities of every position strictly
        # between them, accumulated per gap in O(K^2) total.
        for i in range(len(inner)):
            gap_product = 1.0
            for j in range(i + 1, len(inner)):
                add(
                    cut_positions[j] - cut_positions[i],
                    cut_probs[i] * cut_probs[j] * gap_product,
                )
                gap_product *= cut_stays[j]

        add(n, no_cut_prob)
        return expected, no_cut_prob

    # Circular molecule.  Cut positions are residues modulo n.
    prefix, suffix = _prefix_suffix_products(stays)
    no_cut_prob = prefix[-1]
    if not positions:
        return {n: 1.0}, 1.0

    # An arc of length L from cut i to the next cut clockwise appears
    # iff both ends cut and no candidate position strictly inside the
    # arc cuts.  Doubling the positions (and their probabilities)
    # covers arcs that wrap the origin; all products are accumulated
    # in O(K^2).
    doubled = positions + [p + n for p in positions]
    doubled_q = probabilities + probabilities
    doubled_stay = stays + stays
    for i in range(len(positions)):
        gap_product = 1.0
        q = probabilities[i]
        for step in range(1, len(positions)):
            j = i + step
            add(doubled[j] - doubled[i], q * doubled_q[j] * gap_product)
            gap_product *= doubled_stay[j]

    # Full-length fragments: the no-cut outcome and, for each position,
    # the outcome where exactly that one position cuts (linearisation).
    exactly_one = sum(
        q * prefix[i] * suffix[i + 1]
        for i, q in enumerate(probabilities)
    )
    add(n, no_cut_prob + exactly_one)
    return expected, no_cut_prob
