"""Gel-panel design: choose the smallest enzyme set telling all samples apart."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

from .digest import digest, enzyme_cut_positions


@dataclass(frozen=True)
class PanelSample:
    """One DNA sample to fingerprint."""

    name: str
    sequence: str
    circular: bool


@dataclass(frozen=True)
class PanelEnzyme:
    """One candidate restriction enzyme."""

    name: str
    site: str
    cut_offset: int


@dataclass(frozen=True)
class PanelResult:
    """Outcome of a panel design.

    ``enzymes`` holds the selected enzyme names (sorted); it is empty when no
    subset of the candidates separates every sample pair, in which case
    ``undistinguishable_pairs`` lists the pairs no candidate tells apart.
    """

    enzymes: list[str]
    lanes: dict[str, dict[str, list[int]]]
    pair_enzymes: list[tuple[tuple[str, str], list[str]]]
    undistinguishable_pairs: list[tuple[str, str]]


def single_digest_bands(
    sequence: str, circular: bool, site: str, cut_offset: int
) -> list[int]:
    """Ascending fragment lengths (duplicate bands kept) of a single-enzyme digest."""
    cuts = enzyme_cut_positions(sequence, site, cut_offset, circular)
    return sorted(f.length for f in digest(sequence, cuts, circular))


def bands_distinguishable(
    bands_a: list[int], bands_b: list[int], tolerance_bp: int
) -> bool:
    """True unless both lanes have the same band count and every pair of
    corresponding lengths differs by at most ``tolerance_bp``."""
    if len(bands_a) != len(bands_b):
        return True
    return any(abs(a - b) > tolerance_bp for a, b in zip(bands_a, bands_b))


def minimal_enzyme_panel(
    distinguishing: dict[str, set[int]], pair_total: int
) -> list[str] | None:
    """Smallest enzyme subset whose distinguished pairs cover all pairs.

    Ties are broken by the lexicographically smallest sorted name list
    (combinations over the sorted names are enumerated in that order).
    Returns the sorted names, or None when even all candidates fall short.
    """
    all_pairs = set(range(pair_total))
    names = sorted(distinguishing)
    for size in range(1, len(names) + 1):
        for combo in combinations(names, size):
            covered: set[int] = set()
            for name in combo:
                covered |= distinguishing[name]
            if covered == all_pairs:
                return list(combo)
    return None


def design_panel(
    samples: list[PanelSample],
    enzymes: list[PanelEnzyme],
    tolerance_bp: int,
) -> PanelResult:
    """Single-digest every sample with every candidate enzyme (each enzyme
    forms its own lane), then pick the smallest panel distinguishing all
    sample pairs."""
    # bands[sample][enzyme] -> ascending fragment lengths, duplicates kept
    bands = {
        sample.name: {
            enzyme.name: single_digest_bands(
                sample.sequence, sample.circular, enzyme.site, enzyme.cut_offset
            )
            for enzyme in enzymes
        }
        for sample in samples
    }

    pairs = list(combinations([s.name for s in samples], 2))
    # distinguishing[enzyme] -> indices of the sample pairs it tells apart
    distinguishing: dict[str, set[int]] = {e.name: set() for e in enzymes}
    for index, (name_a, name_b) in enumerate(pairs):
        for enzyme in enzymes:
            if bands_distinguishable(
                bands[name_a][enzyme.name], bands[name_b][enzyme.name], tolerance_bp
            ):
                distinguishing[enzyme.name].add(index)

    selected = minimal_enzyme_panel(distinguishing, len(pairs))
    if selected is None:
        return PanelResult(
            enzymes=[],
            lanes={},
            pair_enzymes=[],
            undistinguishable_pairs=[
                pair
                for index, pair in enumerate(pairs)
                if not any(index in covered for covered in distinguishing.values())
            ],
        )

    return PanelResult(
        enzymes=selected,
        lanes={
            sample.name: {name: bands[sample.name][name] for name in selected}
            for sample in samples
        },
        pair_enzymes=[
            (pair, [name for name in selected if index in distinguishing[name]])
            for index, pair in enumerate(pairs)
        ],
        undistinguishable_pairs=[],
    )
