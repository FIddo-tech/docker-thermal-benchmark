"""
Finite-difference reference solver for 1D steady non-linear heat conduction:
    d/dx( k(T) dT/dx ) + q(x) = 0,   k(T) = k0 (1 + alpha T),   q = Q0 sin(pi x / L)
"""

import json
import os
import numpy as np
from scipy.optimize import root

OUT_DIR = os.environ.get("BENCH_OUTPUT_DIR", "/app/output")

L, N = 1.0, 201
K0, ALPHA, Q0 = 10.0, 0.01, 5000.0
T_LEFT, T_RIGHT = 300.0, 320.0


def k_func(T):
    return K0 * (1.0 + ALPHA * T)


def solve_heat_conduction():
    x = np.linspace(0.0, L, N)
    dx = x[1] - x[0]
    q_src = Q0 * np.sin(np.pi * x / L)

    def residuals(T):
        res = np.empty_like(T)
        res[0], res[-1] = T[0] - T_LEFT, T[-1] - T_RIGHT
        k_plus = k_func(0.5 * (T[2:] + T[1:-1]))
        k_minus = k_func(0.5 * (T[1:-1] + T[:-2]))
        flux_plus = k_plus * (T[2:] - T[1:-1]) / dx
        flux_minus = k_minus * (T[1:-1] - T[:-2]) / dx
        res[1:-1] = (flux_plus - flux_minus) / dx + q_src[1:-1]
        return res

    sol = root(residuals, np.linspace(T_LEFT, T_RIGHT, N), method="hybr", tol=1e-12)
    if not sol.success:
        raise RuntimeError("Solver failed to converge")
    T = sol.x

    # Second-order one-sided gradients at the boundaries; flux q = -k dT/dx
    dTdx_left = (-3 * T[0] + 4 * T[1] - T[2]) / (2 * dx)
    dTdx_right = (3 * T[-1] - 4 * T[-2] + T[-3]) / (2 * dx)
    q_left = -k_func(T[0]) * dTdx_left
    q_right = -k_func(T[-1]) * dTdx_right

    # Trapezoid rule written out by hand so it works on every NumPy version
    total_generated = float(np.sum(0.5 * (q_src[1:] + q_src[:-1]) * np.diff(x)))
    balance_error = abs(total_generated - (q_right - q_left))

    results = {
        "max_temperature": float(np.max(T)),
        "heat_flux_left": float(q_left),
        "heat_flux_right": float(q_right),
        "energy_conservation_error": float(balance_error),
    }
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(os.path.join(OUT_DIR, "results.json"), "w") as f:
        json.dump(results, f, indent=2)
    print("Solver completed. Results written to", OUT_DIR)


if __name__ == "__main__":
    solve_heat_conduction()
