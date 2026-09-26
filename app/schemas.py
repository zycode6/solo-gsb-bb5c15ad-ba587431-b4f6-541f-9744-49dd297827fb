"""Pydantic request/response models for the digest API."""

from __future__ import annotations

import math
import re
from typing import Literal

from pydantic import BaseModel, Field, field_validator, model_validator

from .digest import IUPAC_MATCH

_SEQUENCE_RE = re.compile(r"^[ACGT]+$")
_SITE_RE = re.compile(r"^[ACGTRYSWKMBDHVN]+$")


class EnzymeRule(BaseModel):
    """One restriction enzyme: unique name, IUPAC recognition string and the
    in-string cut offset k (0 <= k <= len(site))."""

    name: str = Field(min_length=1, examples=["EcoRI"])
    site: str = Field(min_length=1, examples=["GAATTC"])
    cut_offset: int = Field(ge=0, examples=[1])

    @field_validator("site")
    @classmethod
    def site_must_be_iupac(cls, value: str) -> str:
        value = value.upper()
        if not _SITE_RE.fullmatch(value):
            raise ValueError(
                "site must contain only IUPAC codes: " + "".join(sorted(IUPAC_MATCH))
            )
        return value

    @model_validator(mode="after")
    def offset_within_site(self) -> "EnzymeRule":
        if self.cut_offset > len(self.site):
            raise ValueError(
                f"cut_offset ({self.cut_offset}) must be <= len(site) ({len(self.site)})"
            )
        return self


class DigestRequest(BaseModel):
    sequence: str = Field(min_length=1, examples=["ACGTGAATTCGGGAATTCAA"])
    topology: Literal["linear", "circular"] = "linear"
    enzymes: list[EnzymeRule] = Field(default_factory=list)

    @field_validator("sequence")
    @classmethod
    def sequence_must_be_acgt(cls, value: str) -> str:
        value = value.upper()
        if not _SEQUENCE_RE.fullmatch(value):
            raise ValueError("sequence must contain only A, C, G, T")
        return value

    @model_validator(mode="after")
    def check_enzymes(self) -> "DigestRequest":
        names = [e.name for e in self.enzymes]
        if len(names) != len(set(names)):
            raise ValueError("enzyme names must be unique")
        for e in self.enzymes:
            if len(e.site) > len(self.sequence):
                raise ValueError(
                    f"recognition site of enzyme {e.name!r} "
                    f"({len(e.site)} bp) is longer than the sequence "
                    f"({len(self.sequence)} bp)"
                )
        return self


class EnzymeCuts(BaseModel):
    name: str
    cuts: list[int]


class FragmentOut(BaseModel):
    start: int
    length: int
    sequence: str


class DigestResponse(BaseModel):
    cuts: list[EnzymeCuts]
    fragments: list[FragmentOut]


MAX_PARTIAL_SEQUENCE_LENGTH = 2000


class PartialEnzymeRule(EnzymeRule):
    """An enzyme rule with a per-site cut probability."""

    p: float = Field(examples=[0.5])

    @field_validator("p")
    @classmethod
    def p_must_be_a_probability(cls, value: float) -> float:
        # Hand-rolled (rather than ge/le) because JSON NaN/Infinity
        # literals parse to non-finite floats that make the standard
        # constraint error body non-serialisable (HTTP 500).
        if not math.isfinite(value) or not 0.0 <= value <= 1.0:
            raise ValueError("p must be a finite number in [0, 1]")
        return value


class PartialDigestRequest(DigestRequest):
    """/partial-digest input: like /digest, with a per-enzyme cut
    probability and a sequence of at most 2000 bp."""

    sequence: str = Field(min_length=1, max_length=MAX_PARTIAL_SEQUENCE_LENGTH)
    enzymes: list[PartialEnzymeRule] = Field(default_factory=list)


class FragmentExpectation(BaseModel):
    """Expected number of fragments of one length (counts may be fractional)."""

    length: int
    count: float


class PartialDigestResponse(BaseModel):
    # Ascending fragment length; same lengths are summed; zeros omitted.
    expected_fragments: list[FragmentExpectation]
    no_cut_probability: float


class PanelSampleIn(BaseModel):
    """One DNA sample: unique name, ACGT sequence and topology."""

    name: str = Field(min_length=1, examples=["wild-type"])
    sequence: str = Field(min_length=1, examples=["AAGAATTCTT"])
    topology: Literal["linear", "circular"] = "linear"

    @field_validator("sequence")
    @classmethod
    def sequence_must_be_acgt(cls, value: str) -> str:
        value = value.upper()
        if not _SEQUENCE_RE.fullmatch(value):
            raise ValueError("sequence must contain only A, C, G, T")
        return value


class DesignPanelRequest(BaseModel):
    """Panel design input: 2-12 uniquely named samples, 1-12 candidate
    enzymes (same rules as /digest) and a non-negative band tolerance."""

    samples: list[PanelSampleIn] = Field(min_length=2, max_length=12)
    enzymes: list[EnzymeRule] = Field(min_length=1, max_length=12)
    tolerance_bp: int = Field(ge=0, examples=[0])

    @model_validator(mode="after")
    def check_panel(self) -> "DesignPanelRequest":
        names = [s.name for s in self.samples]
        if len(names) != len(set(names)):
            raise ValueError("sample names must be unique")
        enzyme_names = [e.name for e in self.enzymes]
        if len(enzyme_names) != len(set(enzyme_names)):
            raise ValueError("enzyme names must be unique")
        for s in self.samples:
            for e in self.enzymes:
                if len(e.site) > len(s.sequence):
                    raise ValueError(
                        f"recognition site of enzyme {e.name!r} "
                        f"({len(e.site)} bp) is longer than the sequence "
                        f"of sample {s.name!r} ({len(s.sequence)} bp)"
                    )
        return self


class SampleLaneBands(BaseModel):
    """One sample's bands in each selected lane (enzyme)."""

    sample: str
    bands: dict[str, list[int]]


class SamplePair(BaseModel):
    samples: list[str] = Field(min_length=2, max_length=2)


class PairDistinction(SamplePair):
    """The selected enzymes able to tell this sample pair apart."""

    enzymes: list[str]


class DesignPanelResponse(BaseModel):
    enzymes: list[str]
    lanes: list[SampleLaneBands]
    pairs: list[PairDistinction]
    undistinguishable_pairs: list[SamplePair]
