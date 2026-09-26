"""Tests for the probabilistic /partial-digest endpoint and its core math."""

from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from app.digest import cut_probabilities, partial_digest_distribution
from app.main import app

client = TestClient(app)


# ---------- core distribution: hand-computed cases ----------


def test_empty_cut_set_means_whole_molecule():
    for circular in (False, True):
        expected, no_cut = partial_digest_distribution(10, {}, circular)
        assert expected == {10: 1.0}
        assert no_cut == 1.0


def test_linear_single_site_half():
    # One inner cut with p=0.5 on a length-10 linear molecule.
    expected, no_cut = partial_digest_distribution(10, {4: 0.5}, circular=False)
    assert no_cut == pytest.approx(0.5)
    assert expected == {4: 0.5, 6: 0.5, 10: 0.5}


def test_linear_two_sites_independent():
    # Cuts at 4 and 8, p=0.5 each.
    expected, no_cut = partial_digest_distribution(10, {4: 0.5, 8: 0.5}, False)
    assert no_cut == pytest.approx(0.25)
    # lengths 4,8: 0.25 | 4,4,2: 0.25 | 8,2: 0.25 | 10: 0.25
    assert expected == {
        2: pytest.approx(0.5),
        4: pytest.approx(0.75),
        6: pytest.approx(0.25),
        8: pytest.approx(0.25),
        10: pytest.approx(0.25),
    }


def test_linear_ends_are_not_random_cuts():
    # Cut positions 0 and n (strand ends) never count as nicks, for any p.
    for ends in ({0: 0.5, 10: 0.5}, {0: 1.0, 10: 1.0}, {0: 0.0, 10: 0.0}):
        expected, no_cut = partial_digest_distribution(10, ends, circular=False)
        assert expected == {10: 1.0}
        assert no_cut == 1.0


def test_circular_single_site_full_length_either_way():
    # No cut -> one full-length fragment; exactly one cut -> one
    # full-length linearised fragment.  Their expected counts add, the
    # no-cut probability counts only the former.
    expected, no_cut = partial_digest_distribution(10, {3: 0.5}, circular=True)
    assert no_cut == pytest.approx(0.5)
    assert expected == {10: 1.0}


def test_circular_two_sites_half():
    # Two sites opposite-ish at 4 and 8 on a 10-bp circle.
    expected, no_cut = partial_digest_distribution(10, {4: 0.5, 8: 0.5}, True)
    assert no_cut == pytest.approx(0.25)
    # no cut: 0.25 of a length-10 fragment; one cut: 2 * 0.25 = 0.5 more.
    assert expected[10] == pytest.approx(0.75)
    # two cuts: arcs 4 and 6, each with probability 0.25.
    assert expected[4] == pytest.approx(0.25)
    assert expected[6] == pytest.approx(0.25)


def test_cut_probabilities_aggregate_shared_position():
    probs = cut_probabilities([([2, 5], 0.5), ([5], 0.25)])
    assert probs[2] == pytest.approx(0.5)
    assert probs[5] == pytest.approx(1.0 - 0.5 * 0.75)


def test_cut_probabilities_zero_probability_enzyme_ignored():
    assert cut_probabilities([([1], 0.0)]) == {}
    # p=1 dominates regardless of other enzymes at the same position.
    assert cut_probabilities([([2], 1.0), ([2], 0.3)]) == {2: 1.0}


# ---------- HTTP API ----------


def test_partial_digest_linear_hand_computed():
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "AAGAATTCAA",  # length 10
            "topology": "linear",
            "enzymes": [
                {"name": "EcoRI", "site": "GAATTC", "cut_offset": 1, "p": 0.5}
            ],
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    # EcoRI on this site cuts at 3 and 7 -> segments 3, 4, 3.
    assert body["no_cut_probability"] == pytest.approx(0.25)
    by_length = {f["length"]: f["count"] for f in body["expected_fragments"]}
    # length 3 is produced by either single cut (0.25 each) and twice
    # when both cuts fire (2 * 0.25) = 1.0.
    assert by_length[3] == pytest.approx(1.0)
    # length 4 only when both positions cut; length 7 with cut at 3 only.
    assert by_length[4] == pytest.approx(0.25)
    assert by_length[7] == pytest.approx(0.5)
    assert by_length[10] == pytest.approx(0.25)


def test_partial_digest_response_sorted_by_length_and_zeros_omitted():
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "ACGTACGTAC",  # 10
            "topology": "linear",
            "enzymes": [
                {"name": "E", "site": "ACGT", "cut_offset": 2, "p": 0.5}
            ],
        },
    )
    assert resp.status_code == 200
    lengths = [f["length"] for f in resp.json()["expected_fragments"]]
    assert lengths == sorted(lengths)
    assert all(f["count"] > 0 for f in resp.json()["expected_fragments"])


def test_partial_digest_two_enzymes_same_position():
    # Two identical enzymes at the same unique site combine to
    # 1 - (1-0.5)(1-0.5) = 0.75 per-site cut probability.
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "GAATTCAA",
            "enzymes": [
                {"name": "A", "site": "GAATTC", "cut_offset": 1, "p": 0.5},
                {"name": "B", "site": "GAATTC", "cut_offset": 1, "p": 0.5},
            ],
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    # Each enzyme cuts at 1 and 5 independently.
    q = 0.75
    assert body["no_cut_probability"] == pytest.approx((1 - q) ** 2)
    by_length = {f["length"]: f["count"] for f in body["expected_fragments"]}
    # prefix/suffix/interior expectations
    assert by_length[1] == pytest.approx(q * (1 - q) + q * q)
    assert by_length[4] == pytest.approx(q * q)
    assert by_length[7] == pytest.approx(q * (1 - q))
    assert by_length[8] == pytest.approx((1 - q) ** 2)


def test_partial_digest_empty_enzyme_list_does_not_cut():
    resp = client.post(
        "/partial-digest",
        json={"sequence": "ACGT", "topology": "linear", "enzymes": []},
    )
    assert resp.status_code == 200
    assert resp.json() == {
        "expected_fragments": [{"length": 4, "count": 1.0}],
        "no_cut_probability": 1.0,
    }


def test_partial_digest_p_zero_enzyme_does_not_cut():
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "GAATTC",
            "topology": "circular",
            "enzymes": [
                {"name": "E", "site": "GAATTC", "cut_offset": 1, "p": 0.0}
            ],
        },
    )
    assert resp.status_code == 200
    assert resp.json() == {
        "expected_fragments": [{"length": 6, "count": 1.0}],
        "no_cut_probability": 1.0,
    }


def test_partial_digest_p_one_matches_deterministic_digest():
    # p=1 circular with two cuts gives exactly the /digest fragments.
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "GAATTC",
            "topology": "circular",
            "enzymes": [
                {"name": "EcoRI", "site": "GAATTC", "cut_offset": 1, "p": 1.0}
            ],
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["no_cut_probability"] == 0.0
    assert body["expected_fragments"] == [
        {"length": 2, "count": 1.0},
        {"length": 4, "count": 1.0},
    ]


def test_partial_digest_p_one_linear_with_end_cut_has_no_end_nick():
    # EcoRI cuts at 1 and 5 in an 8-bp linear sequence; p=1 everywhere.
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "GAATTCAA",
            "enzymes": [
                {"name": "EcoRI", "site": "GAATTC", "cut_offset": 1, "p": 1.0}
            ],
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["no_cut_probability"] == 0.0
    assert body["expected_fragments"] == [
        {"length": 1, "count": 1.0},
        {"length": 3, "count": 1.0},
        {"length": 4, "count": 1.0},
    ]


def test_partial_digest_default_topology_linear():
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "acgt",
            "enzymes": [
                {"name": "X", "site": "ac", "cut_offset": 0, "p": 1.0}
            ],
        },
    )
    assert resp.status_code == 200
    # Linear cut at 0 is the molecule end -> no nick at all.
    assert resp.json() == {
        "expected_fragments": [{"length": 4, "count": 1.0}],
        "no_cut_probability": 1.0,
    }


# ---------- validation (HTTP 422) ----------


def test_partial_digest_reject_probability_above_one():
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "ACGT",
            "enzymes": [{"name": "X", "site": "AC", "cut_offset": 0, "p": 1.1}],
        },
    )
    assert resp.status_code == 422


def test_partial_digest_reject_negative_probability():
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "ACGT",
            "enzymes": [{"name": "X", "site": "AC", "cut_offset": 0, "p": -0.01}],
        },
    )
    assert resp.status_code == 422


@pytest.mark.parametrize("token", ["NaN", "Infinity", "-Infinity"])
def test_partial_digest_reject_non_finite_probability(token):
    # Raw JSON, since Python's json module refuses NaN/Infinity by default.
    raw = (
        '{"sequence":"ACGT","enzymes":['
        '{"name":"X","site":"AC","cut_offset":0,'
        f'"p":{token}'
        "}]}"
    )
    resp = client.post(
        "/partial-digest",
        content=raw,
        headers={"Content-Type": "application/json"},
    )
    assert resp.status_code == 422


def test_partial_digest_reject_sequence_over_2000bp():
    resp = client.post(
        "/partial-digest",
        json={"sequence": "A" * 2001, "enzymes": []},
    )
    assert resp.status_code == 422


def test_partial_digest_accepts_exactly_2000bp():
    resp = client.post(
        "/partial-digest",
        json={"sequence": "A" * 2000, "enzymes": []},
    )
    assert resp.status_code == 200
    assert resp.json()["expected_fragments"] == [{"length": 2000, "count": 1.0}]


def test_partial_digest_inherits_base_validation():
    # non-ACGT sequence
    assert (
        client.post(
            "/partial-digest",
            json={
                "sequence": "ACGX",
                "enzymes": [{"name": "X", "site": "AC", "cut_offset": 0, "p": 0.5}],
            },
        ).status_code
        == 422
    )
    # duplicate enzyme names
    assert (
        client.post(
            "/partial-digest",
            json={
                "sequence": "ACGT",
                "enzymes": [
                    {"name": "X", "site": "AC", "cut_offset": 0, "p": 0.5},
                    {"name": "X", "site": "GT", "cut_offset": 1, "p": 0.5},
                ],
            },
        ).status_code
        == 422
    )
    # site longer than sequence
    assert (
        client.post(
            "/partial-digest",
            json={
                "sequence": "ACGT",
                "enzymes": [
                    {"name": "X", "site": "GAATTC", "cut_offset": 1, "p": 0.5}
                ],
            },
        ).status_code
        == 422
    )
