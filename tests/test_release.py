"""Tests for the /release-fragment endpoint and the release-design logic."""

from fastapi.testclient import TestClient

from app.main import app
from app.release import ReleaseEnzyme, design_release, target_fragment

client = TestClient(app)

ECORI = ReleaseEnzyme("EcoRI", "GAATTC", 1)
BAMHI = ReleaseEnzyme("BamHI", "GGATCC", 1)


# ---------- core logic ----------


def test_target_fragment_linear():
    frag = target_fragment("AAGAATTCTT", False, 3, 7)
    assert (frag.start, frag.length, frag.sequence) == (3, 4, "AATT")


def test_target_fragment_circular_wraps_origin():
    frag = target_fragment("GAATTC", True, 5, 1)
    assert (frag.start, frag.length, frag.sequence) == (5, 2, "CG")


def test_single_enzyme_releases_fragment():
    # EcoRI cuts AAGAATTCTTGGATCCAA at 3 and 7; BamHI at 11 and 15.
    result = design_release("AAGAATTCTTGGATCCAA", False, 3, 7, [ECORI, BAMHI])
    assert result is not None
    assert result.enzymes == ["EcoRI"]
    assert result.cuts == [3, 7]


def test_two_enzymes_needed_for_two_ends():
    # BamHI cuts 2/6, EcoRI cuts 8/12: target [6, 8) needs one cut from each.
    result = design_release("AGGATCCGAATTCA", False, 6, 8, [BAMHI, ECORI])
    assert result is not None
    assert result.enzymes == ["BamHI", "EcoRI"]
    assert result.cuts == [2, 6, 8, 12]


def test_interior_cut_disqualifies_enzyme():
    # EcoRI cuts 5/9/13/17; target [5, 17) has interior cuts at 9 and 13.
    assert design_release("ACGTGAATTCGGGAATTCAA", False, 5, 17, [ECORI]) is None


def test_same_position_multiple_enzymes_count_as_one_cut():
    # Both enzymes cut exactly {3, 7}: the combined cut list is de-duplicated
    # and the lexicographically smaller name wins the tie.
    result = design_release(
        "AAGAATTCTT",
        False,
        3,
        7,
        [ReleaseEnzyme("EcoRI-like", "GAATTC", 1), ECORI],
    )
    assert result is not None
    assert result.enzymes == ["EcoRI"]
    assert result.cuts == [3, 7]


def test_lexicographic_tie_break_among_equal_sets():
    # "Z" and "A" both cut exactly {3, 7}; the sorted name list decides.
    result = design_release(
        "AAGAATTCTT",
        False,
        3,
        7,
        [ReleaseEnzyme("Z", "GAATTC", 1), ReleaseEnzyme("A", "AATT", 0)],
    )
    assert result is not None
    assert result.enzymes == ["A"]
    assert result.cuts == [3, 7]


def test_linear_natural_ends_need_no_cut_there():
    # start=0 is a natural end; only the cut at 3 is required.
    result = design_release("AAGAATTCTT", False, 0, 3, [ECORI])
    assert result is not None
    assert result.enzymes == ["EcoRI"]
    assert result.cuts == [3, 7]


def test_full_length_linear_target_needs_no_enzyme():
    # Both ends are natural ends, so the empty subset is minimal even
    # though the candidate enzyme would cut inside.
    result = design_release("AAGAATTCTT", False, 0, 10, [ECORI])
    assert result is not None
    assert result.enzymes == []
    assert result.cuts == []


def test_circular_target_requires_cuts_at_both_ends():
    # Circular GAATTC: EcoRI cuts at 1 and 5; [1, 5) and wrapping [5, 1).
    for start, end in ((1, 5), (5, 1)):
        result = design_release("GAATTC", True, start, end, [ECORI])
        assert result is not None
        assert result.enzymes == ["EcoRI"]
        assert result.cuts == [1, 5]


def test_circular_has_no_natural_ends():
    # Coordinate 0 is not a cut here, so [0, 5) cannot be released.
    assert design_release("GAATTC", True, 0, 5, [ECORI]) is None


def test_circular_interior_checked_along_forward_arc():
    # EcoRI cuts circular GAATTCGAATTC at 1, 5, 7, 11.  The wrapping
    # target [5, 1) contains cuts 7 and 11 -> infeasible; [1, 5) is clean.
    assert design_release("GAATTCGAATTC", True, 5, 1, [ECORI]) is None
    result = design_release("GAATTCGAATTC", True, 1, 5, [ECORI])
    assert result is not None
    assert result.enzymes == ["EcoRI"]


def test_no_candidates_only_full_length_linear_solves():
    assert design_release("AAGAATTCTT", False, 3, 7, []) is None
    assert design_release("GAATTC", True, 1, 5, []) is None
    result = design_release("AAGAATTCTT", False, 0, 10, [])
    assert result is not None
    assert result.enzymes == []


# ---------- HTTP API ----------


def test_release_fragment_single_enzyme():
    resp = client.post(
        "/release-fragment",
        json={
            "sequence": "AAGAATTCTTGGATCCAA",
            "topology": "linear",
            "start": 3,
            "end": 7,
            "enzymes": [
                {"name": "EcoRI", "site": "GAATTC", "cut_offset": 1},
                {"name": "BamHI", "site": "GGATCC", "cut_offset": 1},
            ],
        },
    )
    assert resp.status_code == 200
    assert resp.json() == {
        "solved": True,
        "enzymes": ["EcoRI"],
        "fragment": {"start": 3, "length": 4, "sequence": "AATT"},
        "cuts": [3, 7],
    }


def test_release_fragment_two_enzymes():
    resp = client.post(
        "/release-fragment",
        json={
            "sequence": "AGGATCCGAATTCA",
            "start": 6,
            "end": 8,
            "enzymes": [
                {"name": "BamHI", "site": "GGATCC", "cut_offset": 1},
                {"name": "EcoRI", "site": "GAATTC", "cut_offset": 1},
            ],
        },
    )
    assert resp.status_code == 200
    assert resp.json() == {
        "solved": True,
        "enzymes": ["BamHI", "EcoRI"],
        "fragment": {"start": 6, "length": 2, "sequence": "CG"},
        "cuts": [2, 6, 8, 12],
    }


def test_release_fragment_circular_wrap():
    resp = client.post(
        "/release-fragment",
        json={
            "sequence": "GAATTC",
            "topology": "circular",
            "start": 5,
            "end": 1,
            "enzymes": [{"name": "EcoRI", "site": "GAATTC", "cut_offset": 1}],
        },
    )
    assert resp.status_code == 200
    assert resp.json() == {
        "solved": True,
        "enzymes": ["EcoRI"],
        "fragment": {"start": 5, "length": 2, "sequence": "CG"},
        "cuts": [1, 5],
    }


def test_release_fragment_unsolved_returns_clear_result():
    resp = client.post(
        "/release-fragment",
        json={
            "sequence": "ACGTGAATTCGGGAATTCAA",
            "start": 5,
            "end": 17,
            "enzymes": [{"name": "EcoRI", "site": "GAATTC", "cut_offset": 1}],
        },
    )
    assert resp.status_code == 200
    assert resp.json() == {"solved": False, "enzymes": [], "fragment": None, "cuts": []}


def test_release_fragment_defaults_to_linear_and_lowercase_ok():
    resp = client.post(
        "/release-fragment",
        json={
            "sequence": "aagaattctt",
            "start": 3,
            "end": 7,
            "enzymes": [{"name": "ecori", "site": "gaattc", "cut_offset": 1}],
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["solved"] is True
    assert body["enzymes"] == ["ecori"]
    assert body["fragment"] == {"start": 3, "length": 4, "sequence": "AATT"}


def test_release_fragment_full_length_linear_target():
    resp = client.post(
        "/release-fragment",
        json={
            "sequence": "AAGAATTCTT",
            "start": 0,
            "end": 10,
            "enzymes": [{"name": "EcoRI", "site": "GAATTC", "cut_offset": 1}],
        },
    )
    assert resp.status_code == 200
    assert resp.json() == {
        "solved": True,
        "enzymes": [],
        "fragment": {"start": 0, "length": 10, "sequence": "AAGAATTCTT"},
        "cuts": [],
    }


def test_release_fragment_no_candidate_enzymes():
    resp = client.post(
        "/release-fragment",
        json={"sequence": "AAGAATTCTT", "start": 3, "end": 7, "enzymes": []},
    )
    assert resp.status_code == 200
    assert resp.json()["solved"] is False


# ---------- validation (HTTP 422) ----------


def _release_payload(**overrides):
    payload = {
        "sequence": "AAGAATTCTT",
        "topology": "linear",
        "start": 3,
        "end": 7,
        "enzymes": [{"name": "EcoRI", "site": "GAATTC", "cut_offset": 1}],
    }
    payload.update(overrides)
    return payload


def test_reject_bad_linear_interval():
    bad = [(7, 3), (3, 3), (0, 11), (-1, 7)]
    for start, end in bad:
        resp = client.post(
            "/release-fragment", json=_release_payload(start=start, end=end)
        )
        assert resp.status_code == 422, (start, end)


def test_reject_bad_circular_interval():
    # start == end is not a valid circular target.
    resp = client.post(
        "/release-fragment",
        json=_release_payload(topology="circular", start=3, end=3),
    )
    assert resp.status_code == 422
    # Coordinates must be < sequence length (no natural end n on a circle).
    resp = client.post(
        "/release-fragment",
        json=_release_payload(topology="circular", start=3, end=10),
    )
    assert resp.status_code == 422


def test_reject_too_many_enzymes():
    enzymes = [{"name": f"E{i}", "site": "GAATTC", "cut_offset": 1} for i in range(13)]
    resp = client.post("/release-fragment", json=_release_payload(enzymes=enzymes))
    assert resp.status_code == 422


def test_reject_duplicate_enzyme_names():
    enzymes = [
        {"name": "X", "site": "GAATTC", "cut_offset": 1},
        {"name": "X", "site": "GGATCC", "cut_offset": 1},
    ]
    resp = client.post("/release-fragment", json=_release_payload(enzymes=enzymes))
    assert resp.status_code == 422


def test_reject_site_longer_than_sequence():
    enzymes = [{"name": "X", "site": "GAATTCGAATTCG", "cut_offset": 1}]
    resp = client.post("/release-fragment", json=_release_payload(enzymes=enzymes))
    assert resp.status_code == 422


def test_reject_non_acgt_sequence_and_non_iupac_site():
    resp = client.post(
        "/release-fragment", json=_release_payload(sequence="AAGAATTCTX")
    )
    assert resp.status_code == 422
    bad_site = [{"name": "X", "site": "ACZ", "cut_offset": 0}]
    resp = client.post("/release-fragment", json=_release_payload(enzymes=bad_site))
    assert resp.status_code == 422
