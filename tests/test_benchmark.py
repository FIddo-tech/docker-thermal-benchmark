"""
Verifier for the thermal-conduction benchmark.
Checks the solver's output against an EXACT analytic solution (Kirchhoff transform),
not just against numbers the solver reports about itself.
"""

import json
import os
import numpy as np
import pytest

RESULTS_PATH = os.path.join(os.environ.get("BENCH_OUTPUT_DIR", "/app/output"), "results.json")

L, K0, ALPHA, Q0 = 1.0, 10.0, 0.01, 5000.0
T_LEFT, T_RIGHT = 300.0, 320.0


def analytic():
    """With theta = k0 (T + alpha T^2 / 2) the PDE becomes theta'' = -q (linear)."""
    th = lambda T: K0 * (T + ALPHA * T ** 2 / 2)
    th0, thL = th(T_LEFT), th(T_RIGHT)
    x = np.linspace(0, L, 200001)
    theta = th0 + (thL - th0) * x / L + Q0 * L ** 2 / np.pi ** 2 * np.sin(np.pi * x / L)
    T = (-1 + np.sqrt(1 + 2 * ALPHA * theta / K0)) / ALPHA
    slope = (thL - th0) / L
    return {
        "max_temperature": float(T.max()),
        "heat_flux_left": float(-(slope + Q0 * L / np.pi)),
        "heat_flux_right": float(-(slope - Q0 * L / np.pi)),
        "total_generated": 2 * Q0 * L / np.pi,
    }


@pytest.fixture
def results():
    assert os.path.exists(RESULTS_PATH), f"{RESULTS_PATH} does not exist"
    with open(RESULTS_PATH) as f:
        return json.load(f)


def test_json_schema(results):
    for key in ("max_temperature", "heat_flux_left", "heat_flux_right",
                "energy_conservation_error"):
        assert key in results, f"Missing key: {key}"
        assert np.isfinite(results[key])


def test_max_temperature_matches_analytic(results):
    ref = analytic()["max_temperature"]
    assert results["max_temperature"] > T_RIGHT          # genuine interior peak
    assert results["max_temperature"] == pytest.approx(ref, abs=0.05)


def test_boundary_fluxes_match_analytic(results):
    ref = analytic()
    assert results["heat_flux_left"] == pytest.approx(ref["heat_flux_left"], rel=2e-3)
    assert results["heat_flux_right"] == pytest.approx(ref["heat_flux_right"], rel=2e-3)


def test_global_energy_balance_independent(results):
    """Net outward flux must equal the exact integral of q(x) = 2 Q0 L / pi."""
    net_out = results["heat_flux_right"] - results["heat_flux_left"]
    assert net_out == pytest.approx(analytic()["total_generated"], abs=2.0)


def test_reported_energy_error_is_small(results):
    assert results["energy_conservation_error"] < 2.0
