"""Checks the diminishing returns term against Equation 2 of Roper, Ng & Pettit (2023, JTLU).

For a count category, the j-th nearest destination should receive
    U_j = exp(-lambda * (j-1)) - exp(-lambda * j)
so that n destinations at zero distance score 1 - exp(-lambda * n).
These tests need no network or OSM download: they call the scoring functions directly.
"""
import numpy as np
import pandas as pd
import pytest

from there_access.index import _calculate_accessibility_score
from there_access.index_oop import AccessibilityCalculator

LAMBDA = 0.3      # default diminishing returns constant for count categories
GAMMA = 0.001     # default distance decay constant (per metre)
NUM_POIS = 25
CAP = 10_000      # search distance; pandana returns this value when fewer POIs are found


def make_inputs(distances, weights=None):
    """One origin row with the given sorted distances; unused slots padded like pandana does.

    weights gives each destination's attractiveness; None means a count category.
    """
    n = len(distances)
    weights = [1] * n if weights is None else list(weights)
    dist_cols = list(distances) + [CAP] * (NUM_POIS - n)
    id_cols = list(range(n)) + [np.nan] * (NUM_POIS - n)
    access = pd.DataFrame(
        [dist_cols + id_cols],
        columns=list(range(1, NUM_POIS + 1)) + [f"poi{i}" for i in range(1, NUM_POIS + 1)],
    )
    relevant_pois = pd.DataFrame({"attract": weights, "chain_weight": [1] * n}, index=range(n))
    return access, relevant_pois


def equation_2(distances, weights=None, lam=LAMBDA, gamma=GAMMA):
    weights = np.ones(len(distances)) if weights is None else np.asarray(weights, dtype=float)
    cumulative = np.cumsum(weights)
    marginal = np.exp(-lam * (cumulative - weights)) - np.exp(-lam * cumulative)
    return float((marginal * np.exp(-gamma * np.asarray(distances, dtype=float))).sum())


SCORERS = {
    "index": lambda a, p: _calculate_accessibility_score(a, p, GAMMA, LAMBDA, NUM_POIS, CAP),
    "index_oop": lambda a, p: AccessibilityCalculator.calculate_score(a, p, GAMMA, LAMBDA, NUM_POIS, CAP),
}


@pytest.mark.parametrize("scorer", SCORERS)
@pytest.mark.parametrize("n", [1, 3, 25])
def test_destinations_at_zero_distance(scorer, n):
    access, pois = make_inputs([0] * n)
    score = float(SCORERS[scorer](access, pois)[0])
    assert score == pytest.approx(1 - np.exp(-LAMBDA * n))


@pytest.mark.parametrize("scorer", SCORERS)
def test_paper_example_three_destinations(scorer):
    # Section 2.2: three destinations, lambda = 0.3, give a maximum score of about 60%
    access, pois = make_inputs([0, 0, 0])
    assert float(SCORERS[scorer](access, pois)[0]) == pytest.approx(0.593, abs=1e-3)


@pytest.mark.parametrize("scorer", SCORERS)
def test_mixed_distances_match_equation_2(scorer):
    distances = [50, 120, 400, 800, 1250, 2000]
    access, pois = make_inputs(distances)
    assert float(SCORERS[scorer](access, pois)[0]) == pytest.approx(equation_2(distances))


@pytest.mark.parametrize("scorer", SCORERS)
def test_category_can_approach_full_score(scorer):
    # previously capped at exp(-lambda) = 0.741
    access, pois = make_inputs([0] * NUM_POIS)
    assert float(SCORERS[scorer](access, pois)[0]) > 0.99


@pytest.mark.parametrize("scorer", ["index", "index_oop"])
def test_weighted_category_matches_equation_2(scorer):
    # weighted category (eg. job counts): attractiveness enters the diminishing returns term
    lam = 0.002
    distances = [30, 150, 420, 600, 1100, 1900]
    jobs = [120, 35.5, 800, 4, 260.25, 1500]
    access, pois = make_inputs(distances, jobs)
    if scorer == "index":
        score = _calculate_accessibility_score(access, pois, GAMMA, lam, NUM_POIS, CAP)
    else:
        score = AccessibilityCalculator.calculate_score(access, pois, GAMMA, lam, NUM_POIS, CAP)
    assert float(score[0]) == pytest.approx(equation_2(distances, jobs, lam=lam))
