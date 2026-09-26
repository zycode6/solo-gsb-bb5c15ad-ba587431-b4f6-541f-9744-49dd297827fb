"""FastAPI application exposing the restriction-digest simulation."""

from __future__ import annotations

import math

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse

from .digest import (
    cut_probabilities,
    digest,
    enzyme_cut_positions,
    partial_digest_distribution,
)
from .panel import PanelEnzyme, PanelSample, design_panel
from .schemas import (
    DesignPanelRequest,
    DesignPanelResponse,
    DigestRequest,
    DigestResponse,
    EnzymeCuts,
    FragmentExpectation,
    FragmentOut,
    PairDistinction,
    PartialDigestRequest,
    PartialDigestResponse,
    SampleLaneBands,
    SamplePair,
)

app = FastAPI(
    title="DNA Restriction Digest Simulator",
    version="1.0.0",
    description=(
        "Simulates restriction-enzyme digestion of linear or circular DNA. "
        "Recognition sites (IUPAC) and their reverse complements are searched "
        "with overlaps kept; all enzymes cut simultaneously. "
        "Also designs minimal enzyme panels that distinguish DNA samples."
    ),
)


def _json_safe(value: object) -> object:
    """Coerce request-validation error data into JSON-serialisable form:
    non-finite floats (NaN/Infinity from request bodies) are stringified
    and arbitrary context objects (e.g. the original ValueError) dropped."""
    if isinstance(value, float):
        if math.isnan(value):
            return "NaN"
        if math.isinf(value):
            return "Infinity" if value > 0 else "-Infinity"
        return value
    if isinstance(value, dict):
        return {k: _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    return str(value)


@app.exception_handler(RequestValidationError)
async def validation_exception_handler(
    request: Request, exc: RequestValidationError
) -> JSONResponse:
    return JSONResponse(
        status_code=422,
        content={"detail": _json_safe(exc.errors())},
    )


@app.get("/health", tags=["meta"])
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.post("/digest", response_model=DigestResponse, tags=["digest"])
def run_digest(request: DigestRequest) -> DigestResponse:
    circular = request.topology == "circular"

    per_enzyme: list[EnzymeCuts] = []
    all_cuts: list[int] = []
    for enzyme in request.enzymes:
        positions = enzyme_cut_positions(
            request.sequence, enzyme.site, enzyme.cut_offset, circular
        )
        per_enzyme.append(EnzymeCuts(name=enzyme.name, cuts=positions))
        all_cuts.extend(positions)

    fragments = digest(request.sequence, all_cuts, circular)

    return DigestResponse(
        cuts=per_enzyme,
        fragments=[FragmentOut(start=f.start, length=f.length, sequence=f.sequence) for f in fragments],
    )


@app.post(
    "/partial-digest",
    response_model=PartialDigestResponse,
    tags=["digest"],
)
def run_partial_digest(request: PartialDigestRequest) -> PartialDigestResponse:
    circular = request.topology == "circular"

    positions_by_enzyme: list[tuple[list[int], float]] = []
    for enzyme in request.enzymes:
        positions = enzyme_cut_positions(
            request.sequence, enzyme.site, enzyme.cut_offset, circular
        )
        positions_by_enzyme.append((positions, enzyme.p))

    cuts = cut_probabilities(positions_by_enzyme)
    expected_by_length, no_cut_prob = partial_digest_distribution(
        len(request.sequence), cuts, circular
    )

    return PartialDigestResponse(
        expected_fragments=[
            FragmentExpectation(length=length, count=expected_by_length[length])
            for length in sorted(expected_by_length)
        ],
        no_cut_probability=no_cut_prob,
    )


@app.post("/design-panel", response_model=DesignPanelResponse, tags=["design-panel"])
def run_design_panel(request: DesignPanelRequest) -> DesignPanelResponse:
    result = design_panel(
        samples=[
            PanelSample(
                name=s.name,
                sequence=s.sequence,
                circular=s.topology == "circular",
            )
            for s in request.samples
        ],
        enzymes=[
            PanelEnzyme(name=e.name, site=e.site, cut_offset=e.cut_offset)
            for e in request.enzymes
        ],
        tolerance_bp=request.tolerance_bp,
    )

    return DesignPanelResponse(
        enzymes=result.enzymes,
        lanes=[
            SampleLaneBands(sample=name, bands=bands)
            for name, bands in result.lanes.items()
        ],
        pairs=[
            PairDistinction(samples=list(pair), enzymes=enzymes)
            for pair, enzymes in result.pair_enzymes
        ],
        undistinguishable_pairs=[
            SamplePair(samples=list(pair)) for pair in result.undistinguishable_pairs
        ],
    )
