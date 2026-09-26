"""Release-fragment design: choose the smallest enzyme set that excises a
target fragment from a linear or circular DNA sequence."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

from .digest import Fragment, enzyme_cut_positions


@dataclass(frozen=True)
class ReleaseEnzyme:
    """One candidate restriction enzyme."""

    name: str
    site: str
    cut_offset: int


@dataclass(frozen=True)
class ReleaseResult:
    """Outcome of a release-fragment design.

    ``enzymes`` holds the selected enzyme names (sorted) and ``cuts`` the
    de-duplicated cut positions they produce together (several enzymes
    cutting at the same position count as one cut point).
    """

    enzymes: list[str]
    cuts: list[int]


def target_fragment(sequence: str, circular: bool, start: int, end: int) -> Fragment:
    """The target fragment: bases from ``start`` to ``end`` along the
    original sequence direction, wrapping the origin on a circle
    (``start != end`` is required there)."""
    if circular:
        n = len(sequence)
        bases = sequence[start:end] if start < end else sequence[start:] + sequence[:end]
        return Fragment(start=start, length=(end - start) % n, sequence=bases)
    return Fragment(start=start, length=end - start, sequence=sequence[start:end])


def _releases_target(
    cuts: set[int], n: int, circular: bool, start: int, end: int
) -> bool:
    """True when cutting at exactly ``cuts`` releases the target fragment:
    both ends are cut points (or the natural ends 0 and n of a linear
    molecule) and no cut falls strictly inside the target.  Cuts outside
    the target are unrestricted."""
    if circular:
        if start not in cuts or end not in cuts:
            return False
        span = (end - start) % n  # 1..n-1 because start != end
        return not any(0 < (c - start) % n < span for c in cuts)
    if start > 0 and start not in cuts:
        return False
    if end < n and end not in cuts:
        return False
    return not any(start < c < end for c in cuts)


def design_release(
    sequence: str,
    circular: bool,
    start: int,
    end: int,
    enzymes: list[ReleaseEnzyme],
) -> ReleaseResult | None:
    """Smallest enzyme subset whose combined digest releases exactly the
    target fragment, or None when no subset of the candidates works.

    Subsets are tried by increasing size and, within one size, in
    lexicographic order of the sorted enzyme-name list (combinations over
    the sorted names), so the first feasible subset found is the answer.
    The empty subset is feasible for a full-length linear target (both
    ends are natural ends).
    """
    cuts_by_name = {
        enzyme.name: enzyme_cut_positions(
            sequence, enzyme.site, enzyme.cut_offset, circular
        )
        for enzyme in enzymes
    }
    names = sorted(cuts_by_name)
    for size in range(0, len(names) + 1):
        for combo in combinations(names, size):
            cuts: set[int] = set()
            for name in combo:
                cuts.update(cuts_by_name[name])
            if _releases_target(cuts, len(sequence), circular, start, end):
                return ReleaseResult(enzymes=list(combo), cuts=sorted(cuts))
    return None
