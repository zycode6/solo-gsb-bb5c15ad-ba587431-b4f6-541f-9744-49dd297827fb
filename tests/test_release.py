"""Tests for the /release-fragment endpoint and the release-design logic."""

from fastapi.testclient import TestClient

from app.main import app
from app.release import ReleaseEnzyme, plan_release, target_fragment

client = TestClient(app)

SEQ = "ACGTGAATTCGGGAATTCAA"  # EcoRI (GAATTC, k=1) cuts [5, 9, 13, 17]
ECORI = {"name": "EcoRI", "site": "GAATTC", "cut_offset": 1}

# L (GAATTC, k=3) cuts [5]; R (GGATCC, k=3) cuts [13]; D (CAAG, k=2) cuts [9].
SEQ2 = "CCGAATTCAAGGATCCTT"
ENZ_L = {"name": "L", "site": "GAATTC", "cut_offset": 3}
ENZ_R = {"name": "R", "site": "GGATCC", "cut_offset": 3}
ENZ_D = {"name": "D", "site": "CAAG", "cut_offset": 2}

# L (GAATTC, k=3) cuts [10, 30]; R (ATTCTG, k=1) cuts [20, 30]: the cut at
# 30 is shared by both enzymes.
SEQ_SHARED = "AAAAAAAGAATTCCCCCCCATTCTGCAGAATTCAAA"
ENZ_SL = {"name": "L", "site": "GAATTC", "cut_offset": 3}
ENZ_SR = {"name": "R", "site": "ATTCTG", "cut_offset": 1}


# ---------- core logic ----------


def test_target_fragment_linear_and_circular_wrap():
    frag = target_fragment("ACGT", False, 1, 3)
    assert (frag.start, frag.length, frag.sequence) == (1, 2, "CG")
    frag = target_fragment("ACGT", True, 3, 2)
    assert (frag.start, frag.length, frag.sequence) == (3, 3, "TAC")


def test_single_enzyme_cutting_both_ends_wins():
    result = plan_release(SEQ, False, 5, 9, [ReleaseEnzyme("EcoRI", "GAATTC", 1)])
    assert result.feasible
    assert result.enzymes == ["EcoRI"]
    assert result.cuts == [5, 9, 13, 17]
    assert (result.fragment.start, result.fragment.length, result.fragment.sequence) == (
        5,
        4,
        "AATT",
    )


def test_enzyme_cutting_inside_target_is_unusable():
    # D cuts at 9 (inside [5, 13)) so it can never be selected; L and R each
    # provide one end, so both are needed.
    enzymes = [ReleaseEnzyme(**e) for e in (ENZ_L, ENZ_R, ENZ_D)]
    result = plan_release(SEQ2, False, 5, 13, enzymes)
    assert result.feasible
    assert result.enzymes == ["L", "R"]
    assert result.cuts == [5, 13]
    assert result.fragment.sequence == "TTCAAGGA"


def test_lexicographic_tie_break():
    enzymes = [
        ReleaseEnzyme("Zeta", "GAATTC", 1),
        ReleaseEnzyme("Alpha", "GAATTC", 1),
    ]
    result = plan_release(SEQ, False, 5, 9, enzymes)
    assert result.enzymes == ["Alpha"]


def test_linear_natural_ends_need_no_cut():
    # start == 0: only the cut at 5 is required.
    result = plan_release(SEQ, False, 0, 5, [ReleaseEnzyme("EcoRI", "GAATTC", 1)])
    assert result.feasible
    assert result.fragment.sequence == "ACGTG"
    # end == n: only the cut at 17 is required.
    result = plan_release(SEQ, False, 17, 20, [ReleaseEnzyme("EcoRI", "GAATTC", 1)])
    assert result.feasible
    assert result.fragment.sequence == "CAA"


def test_full_length_linear_target_needs_no_enzyme():
    enzymes = [
        ReleaseEnzyme("EcoRI", "GAATTC", 1),  # cuts inside -> unusable
        ReleaseEnzyme("None1", "AAAAAA", 1),  # no site -> usable but unneeded
    ]
    result = plan_release(SEQ, False, 0, 20, enzymes)
    assert result.feasible
    assert result.enzymes == []
    assert result.cuts == []
    assert result.fragment.sequence == SEQ


def test_same_point_many_enzymes_counts_as_one_cut():
    # A1 and A2 cut only at 5 (same physical point), B cuts at 13: A1+A2
    # together still provide just one cut point, so B is required.
    enzymes = [
        ReleaseEnzyme("A1", "GAATTC", 3),
        ReleaseEnzyme("A2", "GAATTC", 3),
        ReleaseEnzyme("B", "GGATCC", 3),
    ]
    result = plan_release(SEQ2, False, 5, 13, enzymes)
    assert result.feasible
    assert result.enzymes == ["A1", "B"]
    assert result.cuts == [5, 13]


def test_shared_cut_point_listed_once():
    enzymes = [ReleaseEnzyme("R", "ATTCTG", 1), ReleaseEnzyme("L", "GAATTC", 3)]
    result = plan_release(SEQ_SHARED, False, 10, 20, enzymes)
    assert result.feasible
    assert result.enzymes == ["L", "R"]
    assert result.cuts == [10, 20, 30]  # 30 is cut by both, reported once
    assert (result.fragment.start, result.fragment.length) == (10, 10)


def test_circular_wrap_fragment():
    result = plan_release("GAATTC", True, 5, 1, [ReleaseEnzyme("EcoRI", "GAATTC", 1)])
    assert result.feasible
    assert result.cuts == [1, 5]
    assert (result.fragment.start, result.fragment.length, result.fragment.sequence) == (
        5,
        2,
        "CG",
    )


def test_circular_has_no_natural_ends():
    # The only enzyme cuts at 3 alone; end 1 can never be a cut point.
    result = plan_release("AACGT", True, 3, 1, [ReleaseEnzyme("X", "CGT", 1)])
    assert not result.feasible
    assert result.enzymes == [] and result.fragment is None and result.cuts == []


def test_no_solution_when_ends_not_cuttable():
    result = plan_release(SEQ, False, 6, 8, [ReleaseEnzyme("EcoRI", "GAATTC", 1)])
    assert not result.feasible


def test_no_solution_when_only_candidate_cuts_inside():
    result = plan_release(SEQ, False, 2, 18, [ReleaseEnzyme("EcoRI", "GAATTC", 1)])
    assert not result.feasible


# ---------- HTTP API ----------


def test_release_linear_single_enzyme():
    resp = client.post(
        "/release-fragment",
        json={
            "sequence": SEQ,
            "topology": "linear",
            "start": 5,
            "end": 9,
            "enzymes": [ECORI],
        },
    )
    assert resp.status_code == 200
    assert resp.json() == {
        "feasible": True,
        "enzymes": ["EcoRI"],
        "fragment": {"start": 5, "length": 4, "sequence": "AATT"},
        "cuts": [5, 9, 13, 17],
    }


def test_release_linear_two_enzymes_with_interior_cut_decoy():
    resp = client.post(
        "/release-fragment",
        json={
            "sequence": SEQ2,
            "topology": "linear",
            "start": 5,
            "end": 13,
            "enzymes": [ENZ_L, ENZ_R, ENZ_D],
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["feasible"] is True
    assert body["enzymes"] == ["L", "R"]
    assert body["fragment"] == {"start": 5, "length": 8, "sequence": "TTCAAGGA"}
    assert body["cuts"] == [5, 13]


def test_release_lexicographic_tie_break():
    resp = client.post(
        "/release-fragment",
        json={
            "sequence": SEQ,
            "start": 5,
            "end": 9,
            "enzymes": [
                {"name": "Zeta", "site": "GAATTC", "cut_offset": 1},
                {"name": "Alpha", "site": "GAATTC", "cut_offset": 1},
            ],
        },
    )
    assert resp.status_code == 200
    assert resp.json()["enzymes"] == ["Alpha"]


def test_release_linear_natural_end():
    resp = client.post(
        "/release-fragment",
        json={"sequence": SEQ, "start": 0, "end": 5, "enzymes": [ECORI]},
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["feasible"] is True
    assert body["enzymes"] == ["EcoRI"]
    assert body["fragment"] == {"start": 0, "length": 5, "sequence": "ACGTG"}


def test_release_full_length_linear_target_needs_no_enzyme():
    resp = client.post(
        "/release-fragment",
        json={
            "sequence": SEQ,
            "start": 0,
            "end": 20,
            "enzymes": [ECORI, {"name": "None1", "site": "AAAAAA", "cut_offset": 1}],
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["feasible"] is True
    assert body["enzymes"] == []
    assert body["cuts"] == []
    assert body["fragment"] == {"start": 0, "length": 20, "sequence": SEQ}


def test_release_same_point_counts_once():
    resp = client.post(
        "/release-fragment",
        json={
            "sequence": SEQ2,
            "start": 5,
            "end": 13,
            "enzymes": [
                {"name": "A1", "site": "GAATTC", "cut_offset": 3},
                {"name": "A2", "site": "GAATTC", "cut_offset": 3},
                {"name": "B", "site": "GGATCC", "cut_offset": 3},
            ],
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["enzymes"] == ["A1", "B"]
    assert body["cuts"] == [5, 13]


def test_release_shared_cut_point_deduplicated():
    resp = client.post(
        "/release-fragment",
        json={
            "sequence": SEQ_SHARED,
            "start": 10,
            "end": 20,
            "enzymes": [ENZ_SR, ENZ_SL],
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["enzymes"] == ["L", "R"]
    assert body["cuts"] == [10, 20, 30]
    assert body["fragment"] == {"start": 10, "length": 10, "sequence": "TTCCCCCCCA"}


def test_release_circular_wrap_fragment():
    resp = client.post(
        "/release-fragment",
        json={
            "sequence": "GAATTC",
            "topology": "circular",
            "start": 5,
            "end": 1,
            "enzymes": [ECORI],
        },
    )
    assert resp.status_code == 200
    assert resp.json() == {
        "feasible": True,
        "enzymes": ["EcoRI"],
        "fragment": {"start": 5, "length": 2, "sequence": "CG"},
        "cuts": [1, 5],
    }


def test_release_circular_two_enzymes_wrap():
    resp = client.post(
        "/release-fragment",
        json={
            "sequence": "GGATCCAAAAAAGAATTCAA",
            "topology": "circular",
            "start": 15,
            "end": 3,
            "enzymes": [ENZ_L, ENZ_R],
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["enzymes"] == ["L", "R"]
    assert body["cuts"] == [3, 15]
    assert body["fragment"] == {"start": 15, "length": 8, "sequence": "TTCAAGGA"}


def test_release_circular_end_without_cut_is_infeasible():
    resp = client.post(
        "/release-fragment",
        json={
            "sequence": "AACGT",
            "topology": "circular",
            "start": 3,
            "end": 1,
            "enzymes": [{"name": "X", "site": "CGT", "cut_offset": 1}],
        },
    )
    assert resp.status_code == 200
    assert resp.json() == {"feasible": False, "enzymes": [], "fragment": None, "cuts": []}


def test_release_no_solution():
    resp = client.post(
        "/release-fragment",
        json={"sequence": SEQ, "start": 6, "end": 8, "enzymes": [ECORI]},
    )
    assert resp.status_code == 200
    assert resp.json() == {"feasible": False, "enzymes": [], "fragment": None, "cuts": []}


def test_release_defaults_to_linear_and_lowercase_ok():
    resp = client.post(
        "/release-fragment",
        json={
            "sequence": "acgtgaattcgggaattcaa",
            "start": 5,
            "end": 9,
            "enzymes": [{"name": "EcoRI", "site": "gaattc", "cut_offset": 1}],
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["feasible"] is True
    assert body["fragment"] == {"start": 5, "length": 4, "sequence": "AATT"}


# ---------- validation (HTTP 422) ----------


def _payload(**overrides):
    payload = {
        "sequence": SEQ,
        "topology": "linear",
        "start": 5,
        "end": 9,
        "enzymes": [ECORI],
    }
    payload.update(overrides)
    return payload


def test_reject_linear_bad_interval():
    assert client.post("/release-fragment", json=_payload(start=5, end=5)).status_code == 422
    assert client.post("/release-fragment", json=_payload(start=9, end=5)).status_code == 422
    assert client.post("/release-fragment", json=_payload(start=0, end=21)).status_code == 422
    assert client.post("/release-fragment", json=_payload(start=-1, end=5)).status_code == 422


def test_reject_circular_bad_interval():
    circular = {"topology": "circular"}
    # start == end is not a valid circular target.
    assert (
        client.post("/release-fragment", json=_payload(**circular, start=3, end=3)).status_code
        == 422
    )
    # Circular coordinates must be < sequence length.
    assert (
        client.post("/release-fragment", json=_payload(**circular, start=0, end=20)).status_code
        == 422
    )
    assert (
        client.post("/release-fragment", json=_payload(**circular, start=20, end=3)).status_code
        == 422
    )


def test_reject_missing_or_extra_enzymes():
    assert client.post("/release-fragment", json=_payload(enzymes=[])).status_code == 422
    enzymes = [{"name": f"E{i}", "site": "GAATTC", "cut_offset": 1} for i in range(13)]
    assert client.post("/release-fragment", json=_payload(enzymes=enzymes)).status_code == 422


def test_reject_duplicate_enzyme_names():
    resp = client.post(
        "/release-fragment",
        json=_payload(
            enzymes=[
                {"name": "X", "site": "GAATTC", "cut_offset": 1},
                {"name": "X", "site": "GGATCC", "cut_offset": 1},
            ]
        ),
    )
    assert resp.status_code == 422


def test_reject_non_acgt_sequence():
    resp = client.post("/release-fragment", json=_payload(sequence="ACGTX" * 4))
    assert resp.status_code == 422


def test_reject_non_iupac_site_and_bad_offset():
    bad_site = [{"name": "X", "site": "ACZ", "cut_offset": 0}]
    assert client.post("/release-fragment", json=_payload(enzymes=bad_site)).status_code == 422
    bad_offset = [{"name": "X", "site": "ACG", "cut_offset": 4}]
    assert client.post("/release-fragment", json=_payload(enzymes=bad_offset)).status_code == 422


def test_reject_site_longer_than_sequence():
    resp = client.post(
        "/release-fragment",
        json=_payload(
            sequence="ACGT",
            start=1,
            end=3,
            enzymes=[{"name": "X", "site": "GAATTC", "cut_offset": 1}],
        ),
    )
    assert resp.status_code == 422
