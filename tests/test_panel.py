"""Tests for the /design-panel endpoint and the panel-design logic."""

from fastapi.testclient import TestClient

from app.main import app
from app.panel import (
    bands_distinguishable,
    minimal_enzyme_panel,
    single_digest_bands,
)

client = TestClient(app)


# ---------- core logic ----------


def test_single_digest_bands_are_sorted_and_keep_duplicates():
    bands = single_digest_bands("ACGTGAATTCGGGAATTCAA", False, "GAATTC", 1)
    assert bands == [3, 4, 4, 4, 5]


def test_single_digest_bands_no_cuts_returns_whole_length():
    assert single_digest_bands("ACGT", False, "GG", 0) == [4]
    assert single_digest_bands("ACGT", True, "GG", 0) == [4]


def test_bands_distinguishable_respects_tolerance():
    assert not bands_distinguishable([3, 4, 4], [3, 4, 5], tolerance_bp=1)
    assert bands_distinguishable([3, 4, 4], [3, 4, 5], tolerance_bp=0)
    # Different band counts are always distinguishable.
    assert bands_distinguishable([3, 4], [3, 4, 5], tolerance_bp=10)


def test_minimal_panel_prefers_smallest_set():
    # "C" alone covers both pairs: set size beats lexicographic order.
    panel = minimal_enzyme_panel({"A": {0}, "B": {1}, "C": {0, 1}}, 2)
    assert panel == ["C"]


def test_minimal_panel_lexicographic_tie_break():
    assert minimal_enzyme_panel({"B": {0}, "A": {0}}, 1) == ["A"]
    # First covering size-2 combination in sorted-name order.
    assert minimal_enzyme_panel({"B": {0}, "A": {0}, "C": {1}}, 2) == ["A", "C"]


def test_minimal_panel_returns_none_when_uncoverable():
    assert minimal_enzyme_panel({"A": {0}}, 2) is None


# ---------- HTTP API ----------

# X alone separates s1 from s2/s3, Z alone separates s3 from s1/s2;
# no single enzyme covers all three pairs, and {X, Z} < {Y, Z}.
THREE_SAMPLE_PANEL = {
    "samples": [
        {"name": "s1", "sequence": "AAGAATTCTT"},
        {"name": "s2", "sequence": "AAGGATCCTT"},
        {"name": "s3", "sequence": "TTGGATCCAA"},
    ],
    "enzymes": [
        {"name": "X", "site": "GAATTC", "cut_offset": 1},
        {"name": "Y", "site": "GGATCC", "cut_offset": 1},
        {"name": "Z", "site": "TTGG", "cut_offset": 1},
    ],
    "tolerance_bp": 0,
}


def test_design_panel_minimal_cover_and_pair_report():
    resp = client.post("/design-panel", json=THREE_SAMPLE_PANEL)
    assert resp.status_code == 200
    body = resp.json()
    assert body["enzymes"] == ["X", "Z"]
    assert body["lanes"] == [
        {"sample": "s1", "bands": {"X": [3, 3, 4], "Z": [10]}},
        {"sample": "s2", "bands": {"X": [10], "Z": [10]}},
        {"sample": "s3", "bands": {"X": [10], "Z": [1, 1, 8]}},
    ]
    assert body["pairs"] == [
        {"samples": ["s1", "s2"], "enzymes": ["X"]},
        {"samples": ["s1", "s3"], "enzymes": ["X", "Z"]},
        {"samples": ["s2", "s3"], "enzymes": ["Z"]},
    ]
    assert body["undistinguishable_pairs"] == []


def test_design_panel_single_enzyme_lexicographic_tie():
    resp = client.post(
        "/design-panel",
        json={
            "samples": [
                {"name": "wild-type", "sequence": "AAGAATTCTT"},
                {"name": "mutant", "sequence": "AAGGATCCTT"},
            ],
            "enzymes": [
                {"name": "EcoRI", "site": "GAATTC", "cut_offset": 1},
                {"name": "BamHI", "site": "GGATCC", "cut_offset": 1},
            ],
            "tolerance_bp": 0,
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    # Both enzymes separate the pair alone; sorted-name order breaks the tie.
    assert body["enzymes"] == ["BamHI"]
    assert body["lanes"] == [
        {"sample": "wild-type", "bands": {"BamHI": [10]}},
        {"sample": "mutant", "bands": {"BamHI": [3, 3, 4]}},
    ]
    assert body["pairs"] == [
        {"samples": ["wild-type", "mutant"], "enzymes": ["BamHI"]}
    ]
    assert body["undistinguishable_pairs"] == []


def test_design_panel_tolerance_can_merge_bands():
    payload = {
        "samples": [
            {"name": "s1", "sequence": "AAGAATTCTT"},   # EcoRI bands [3, 3, 4]
            {"name": "s2", "sequence": "AAAGAATTCTT"},  # EcoRI bands [3, 4, 4]
        ],
        "enzymes": [{"name": "EcoRI", "site": "GAATTC", "cut_offset": 1}],
    }
    # A 1 bp difference exceeds tolerance 0 -> distinguishable.
    resp = client.post("/design-panel", json={**payload, "tolerance_bp": 0})
    assert resp.status_code == 200
    assert resp.json()["enzymes"] == ["EcoRI"]

    # Tolerance 1 absorbs the difference -> nothing can tell them apart.
    resp = client.post("/design-panel", json={**payload, "tolerance_bp": 1})
    assert resp.status_code == 200
    body = resp.json()
    assert body["enzymes"] == []
    assert body["lanes"] == []
    assert body["pairs"] == []
    assert body["undistinguishable_pairs"] == [{"samples": ["s1", "s2"]}]


def test_design_panel_reports_only_hopeless_pairs():
    resp = client.post(
        "/design-panel",
        json={
            "samples": [
                {"name": "dup-a", "sequence": "AAGAATTCTT"},
                {"name": "dup-b", "sequence": "AAGAATTCTT"},
                {"name": "other", "sequence": "AAGGATCCTT"},
            ],
            "enzymes": [{"name": "EcoRI", "site": "GAATTC", "cut_offset": 1}],
            "tolerance_bp": 0,
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["enzymes"] == []
    assert body["lanes"] == []
    assert body["pairs"] == []
    assert body["undistinguishable_pairs"] == [{"samples": ["dup-a", "dup-b"]}]


def test_design_panel_honours_per_sample_topology():
    resp = client.post(
        "/design-panel",
        json={
            "samples": [
                {"name": "plasmid", "sequence": "GAATTC", "topology": "circular"},
                {"name": "linear-dna", "sequence": "GAATTC", "topology": "linear"},
            ],
            "enzymes": [{"name": "EcoRI", "site": "GAATTC", "cut_offset": 1}],
            "tolerance_bp": 0,
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["enzymes"] == ["EcoRI"]
    assert body["lanes"] == [
        {"sample": "plasmid", "bands": {"EcoRI": [2, 4]}},
        {"sample": "linear-dna", "bands": {"EcoRI": [1, 1, 4]}},
    ]


def test_design_panel_defaults_to_linear_and_lowercase_ok():
    resp = client.post(
        "/design-panel",
        json={
            "samples": [
                {"name": "s1", "sequence": "aagaattctt"},
                {"name": "s2", "sequence": "aaggatcctt"},
            ],
            "enzymes": [{"name": "EcoRI", "site": "gaattc", "cut_offset": 1}],
            "tolerance_bp": 0,
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["enzymes"] == ["EcoRI"]
    assert body["lanes"][0] == {"sample": "s1", "bands": {"EcoRI": [3, 3, 4]}}


# ---------- validation (HTTP 422) ----------


def _panel_payload(**overrides):
    payload = {
        "samples": [
            {"name": "s1", "sequence": "AAGAATTCTT"},
            {"name": "s2", "sequence": "AAGGATCCTT"},
        ],
        "enzymes": [{"name": "EcoRI", "site": "GAATTC", "cut_offset": 1}],
        "tolerance_bp": 0,
    }
    payload.update(overrides)
    return payload


def test_reject_too_few_samples():
    resp = client.post(
        "/design-panel",
        json=_panel_payload(samples=[{"name": "s1", "sequence": "AAGAATTCTT"}]),
    )
    assert resp.status_code == 422


def test_reject_too_many_samples():
    samples = [{"name": f"s{i}", "sequence": "AAGAATTCTT"} for i in range(13)]
    resp = client.post("/design-panel", json=_panel_payload(samples=samples))
    assert resp.status_code == 422


def test_reject_duplicate_sample_names():
    resp = client.post(
        "/design-panel",
        json=_panel_payload(
            samples=[
                {"name": "s1", "sequence": "AAGAATTCTT"},
                {"name": "s1", "sequence": "AAGGATCCTT"},
            ]
        ),
    )
    assert resp.status_code == 422


def test_reject_missing_or_extra_enzymes():
    assert client.post("/design-panel", json=_panel_payload(enzymes=[])).status_code == 422
    enzymes = [
        {"name": f"E{i}", "site": "GAATTC", "cut_offset": 1} for i in range(13)
    ]
    assert client.post("/design-panel", json=_panel_payload(enzymes=enzymes)).status_code == 422


def test_reject_negative_tolerance():
    resp = client.post("/design-panel", json=_panel_payload(tolerance_bp=-1))
    assert resp.status_code == 422


def test_reject_duplicate_enzyme_names():
    resp = client.post(
        "/design-panel",
        json=_panel_payload(
            enzymes=[
                {"name": "X", "site": "GAATTC", "cut_offset": 1},
                {"name": "X", "site": "GGATCC", "cut_offset": 1},
            ]
        ),
    )
    assert resp.status_code == 422


def test_reject_non_acgt_sample_sequence():
    resp = client.post(
        "/design-panel",
        json=_panel_payload(
            samples=[
                {"name": "s1", "sequence": "AAGAATTCTX"},
                {"name": "s2", "sequence": "AAGGATCCTT"},
            ]
        ),
    )
    assert resp.status_code == 422


def test_reject_non_iupac_site_and_bad_offset():
    bad_site = [{"name": "X", "site": "ACZ", "cut_offset": 0}]
    assert client.post("/design-panel", json=_panel_payload(enzymes=bad_site)).status_code == 422
    bad_offset = [{"name": "X", "site": "ACG", "cut_offset": 4}]
    assert client.post("/design-panel", json=_panel_payload(enzymes=bad_offset)).status_code == 422


def test_reject_site_longer_than_any_sample_sequence():
    resp = client.post(
        "/design-panel",
        json=_panel_payload(
            samples=[
                {"name": "s1", "sequence": "AAGAATTCTT"},
                {"name": "s2", "sequence": "ACGT"},  # shorter than GAATTC (6 bp)
            ]
        ),
    )
    assert resp.status_code == 422


# ---------- /digest regression ----------


def test_digest_endpoint_still_works():
    resp = client.post(
        "/digest",
        json={
            "sequence": "AGGATCCGAATTCA",
            "enzymes": [{"name": "BamHI", "site": "GGATCC", "cut_offset": 1}],
        },
    )
    assert resp.status_code == 200
    assert resp.json()["cuts"] == [{"name": "BamHI", "cuts": [2, 6]}]
