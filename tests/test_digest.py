"""Tests for core digest logic and the HTTP API."""

from fastapi.testclient import TestClient

from app.digest import digest, enzyme_cut_positions, reverse_complement
from app.main import app

client = TestClient(app)


# ---------- core logic ----------


def test_reverse_complement_iupac():
    assert reverse_complement("GAATTC") == "GAATTC"  # palindrome
    assert reverse_complement("AAAGTT") == "AACTTT"
    assert reverse_complement("RGY") == "RCY"


def test_linear_palindromic_enzyme_cuts_both_strands():
    # EcoRI GAATTC, k=1: forward cut s+1, reverse-complement cut s+5.
    cuts = enzyme_cut_positions("ACGTGAATTCGGGAATTCAA", "GAATTC", 1, circular=False)
    assert cuts == [5, 9, 13, 17]


def test_overlapping_matches_kept():
    cuts = enzyme_cut_positions("AAAA", "AA", 1, circular=False)
    assert cuts == [1, 2, 3]


def test_reverse_complement_only_site():
    # Site AAAGTT absent, but its reverse complement AACTTT occurs at index 2.
    cuts = enzyme_cut_positions("TTAACTTTGG", "AAAGTT", 2, circular=False)
    assert cuts == [6]  # 2 + 6 - 2


def test_circular_match_spanning_origin():
    # "TA" wraps around the end of circular ACGT (start 3).
    cuts = enzyme_cut_positions("ACGT", "TA", 1, circular=True)
    assert cuts == [0]  # (3+1) % 4 == (3+2-1) % 4 == 0, de-duplicated


def test_circular_cut_positions_are_modulo_length():
    cuts = enzyme_cut_positions("GAATTC", "GAATTC", 0, circular=True)
    assert cuts == [0]  # (0+0) % 6 == (0+6-0) % 6


def test_linear_fragments():
    frags = digest("ACGTGAATTCGGGAATTCAA", [5, 9, 13, 17], circular=False)
    assert [(f.start, f.length, f.sequence) for f in frags] == [
        (0, 5, "ACGTG"),
        (5, 4, "AATT"),
        (9, 4, "CGGG"),
        (13, 4, "AATT"),
        (17, 3, "CAA"),
    ]


def test_linear_cut_at_boundary_makes_no_empty_fragment():
    frags = digest("ACGT", [0, 4], circular=False)
    assert [(f.start, f.length, f.sequence) for f in frags] == [(0, 4, "ACGT")]


def test_no_cuts_returns_original_sequence():
    for circular in (False, True):
        frags = digest("ACGT", [], circular=circular)
        assert [(f.start, f.length, f.sequence) for f in frags] == [(0, 4, "ACGT")]


def test_circular_fragments_wrap_to_first():
    frags = digest("ACGTACGT", [2, 6], circular=True)
    assert [(f.start, f.length, f.sequence) for f in frags] == [
        (2, 4, "GTAC"),
        (6, 4, "GTAC"),  # wraps: seq[6:] + seq[:2]
    ]


def test_circular_single_cut_linearises_circle():
    frags = digest("GAATTC", [3], circular=True)
    assert [(f.start, f.length, f.sequence) for f in frags] == [(3, 6, "TTCGAA")]


def test_degenerate_iupac_site():
    # BglII RGATCY (R=A/G, Y=C/T) matches AGATCT at index 1.
    cuts = enzyme_cut_positions("CAGATCTG", "RGATCY", 1, circular=False)
    assert cuts == [2, 6]


# ---------- HTTP API ----------


def test_health():
    resp = client.get("/health")
    assert resp.status_code == 200
    assert resp.json() == {"status": "ok"}


def test_digest_endpoint_linear_two_enzymes():
    resp = client.post(
        "/digest",
        json={
            "sequence": "AGGATCCGAATTCA",
            "topology": "linear",
            "enzymes": [
                {"name": "BamHI", "site": "GGATCC", "cut_offset": 1},
                {"name": "EcoRI", "site": "GAATTC", "cut_offset": 1},
            ],
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["cuts"] == [
        {"name": "BamHI", "cuts": [2, 6]},
        {"name": "EcoRI", "cuts": [8, 12]},
    ]
    assert body["fragments"] == [
        {"start": 0, "length": 2, "sequence": "AG"},
        {"start": 2, "length": 4, "sequence": "GATC"},
        {"start": 6, "length": 2, "sequence": "CG"},
        {"start": 8, "length": 4, "sequence": "AATT"},
        {"start": 12, "length": 2, "sequence": "CA"},
    ]


def test_digest_endpoint_circular_wrap_fragment():
    resp = client.post(
        "/digest",
        json={
            "sequence": "GAATTC",
            "topology": "circular",
            "enzymes": [{"name": "EcoRI", "site": "GAATTC", "cut_offset": 1}],
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["cuts"] == [{"name": "EcoRI", "cuts": [1, 5]}]
    assert body["fragments"] == [
        {"start": 1, "length": 4, "sequence": "AATT"},
        {"start": 5, "length": 2, "sequence": "CG"},
    ]


def test_digest_endpoint_defaults_to_linear_and_lowercase_ok():
    resp = client.post(
        "/digest",
        json={
            "sequence": "acgtgaattcaa",
            "enzymes": [{"name": "EcoRI", "site": "gaattc", "cut_offset": 1}],
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["cuts"] == [{"name": "EcoRI", "cuts": [5, 9]}]
    assert body["fragments"][0] == {"start": 0, "length": 5, "sequence": "ACGTG"}


def test_digest_endpoint_no_enzymes_returns_whole_sequence():
    resp = client.post("/digest", json={"sequence": "ACGT", "enzymes": []})
    assert resp.status_code == 200
    assert resp.json() == {
        "cuts": [],
        "fragments": [{"start": 0, "length": 4, "sequence": "ACGT"}],
    }


def test_shared_cut_position_cut_once():
    resp = client.post(
        "/digest",
        json={
            "sequence": "AGAATTCT",
            "enzymes": [
                {"name": "EcoRI", "site": "GAATTC", "cut_offset": 1},
                {"name": "EcoRI-like", "site": "GAATTC", "cut_offset": 1},
            ],
        },
    )
    assert resp.status_code == 200
    frags = resp.json()["fragments"]
    assert [f["start"] for f in frags] == [0, 2, 6]
    assert all(f["length"] > 0 for f in frags)


# ---------- validation (HTTP 422) ----------


def test_reject_duplicate_enzyme_names():
    resp = client.post(
        "/digest",
        json={
            "sequence": "ACGT",
            "enzymes": [
                {"name": "X", "site": "AC", "cut_offset": 0},
                {"name": "X", "site": "GT", "cut_offset": 1},
            ],
        },
    )
    assert resp.status_code == 422


def test_reject_offset_beyond_site_length():
    resp = client.post(
        "/digest",
        json={
            "sequence": "ACGT",
            "enzymes": [{"name": "X", "site": "ACG", "cut_offset": 4}],
        },
    )
    assert resp.status_code == 422


def test_reject_site_longer_than_sequence():
    resp = client.post(
        "/digest",
        json={
            "sequence": "ACGT",
            "enzymes": [{"name": "X", "site": "GAATTC", "cut_offset": 1}],
        },
    )
    assert resp.status_code == 422


def test_reject_non_acgt_sequence():
    resp = client.post(
        "/digest",
        json={
            "sequence": "ACGTX",
            "enzymes": [{"name": "X", "site": "AC", "cut_offset": 0}],
        },
    )
    assert resp.status_code == 422


def test_reject_non_iupac_site():
    resp = client.post(
        "/digest",
        json={
            "sequence": "ACGT",
            "enzymes": [{"name": "X", "site": "ACZ", "cut_offset": 0}],
        },
    )
    assert resp.status_code == 422
