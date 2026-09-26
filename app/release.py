"""Release-fragment design: pick the smallest enzyme set cutting out a target fragment."""

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

    ``feasible`` is False when no subset of the candidate enzymes releases
    exactly the target fragment; then ``enzymes``/``cuts`` are empty and
    ``fragment`` is None.  Otherwise ``enzymes`` holds the selected enzyme
    names (sorted), ``fragment`` the released target fragment and ``cuts``
    the sorted, de-duplicated cut positions the selected enzymes produce
    (several enzymes cutting the same position count as one cut point).
    """

    feasible: bool
    enzymes: list[str]
    fragment: Fragment | None
    cuts: list[int]


def target_fragment(sequence: str, circular: bool, start: int, end: int) -> Fragment:
    """The target fragment: bases from ``start`` to ``end`` in the original
    (5'->3') direction, wrapping around the origin for circular molecules."""
    if circular and end < start:
        bases = sequence[start:] + sequence[:end]
    else:
        bases = sequence[start:end]
    return Fragment(start=start, length=len(bases), sequence=bases)


def plan_release(
    sequence: str,
    circular: bool,
    start: int,
    end: int,
    enzymes: list[ReleaseEnzyme],
) -> ReleaseResult:
    """Smallest enzyme subset whose combined digest releases exactly the
    target fragment.

    Both ends of the target must be cut points of the selected enzymes or
    natural ends of a linear molecule, and no selected enzyme may cut
    strictly inside the target.  Ties are broken by the lexicographically
    smallest sorted name list (combinations over the sorted names are
    enumerated in that order).
    """
    n = len(sequence)
    cuts_by_name = {
        enzyme.name: enzyme_cut_positions(
            sequence, enzyme.site, enzyme.cut_offset, circular
        )
        for enzyme in enzymes
    }

    # Positions strictly inside the target arc: no selected enzyme may cut here.
    if circular:
        arc = (end - start) % n
        interior = {(start + i) % n for i in range(1, arc)}
    else:
        interior = set(range(start + 1, end))
    usable = {
        name: cuts
        for name, cuts in cuts_by_name.items()
        if interior.isdisjoint(cuts)
    }

    # Target ends that must be cut points (linear natural ends need no cut).
    required: set[int] = set()
    if circular or start > 0:
        required.add(start)
    if circular or end < n:
        required.add(end)

    names = sorted(usable)
    for size in range(len(names) + 1):
        for combo in combinations(names, size):
            union: set[int] = set()
            for name in combo:
                union.update(usable[name])
            if required <= union:
                return ReleaseResult(
                    feasible=True,
                    enzymes=list(combo),
                    fragment=target_fragment(sequence, circular, start, end),
                    cuts=sorted(union),
                )
    return ReleaseResult(feasible=False, enzymes=[], fragment=None, cuts=[])
