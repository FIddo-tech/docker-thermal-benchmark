# Task specification

Solve the steady 1D heat conduction equation d/dx( k(T) dT/dx ) + q(x) = 0 on x in [0, L], with k(T) = k0 (1 + alpha T), k0 = 10 W/(m K), alpha = 0.01 1/K, q(x) = Q0 sin(pi x / L), Q0 = 5000 W/m3, L = 1 m. Boundary conditions: T(0) = 300 K and T(L) = 320 K.

The solver must write results.json with:
- max_temperature (K)
- heat_flux_left and heat_flux_right (W/m2, defined as q = -k dT/dx, positive in +x)
- energy_conservation_error (absolute difference between the integral of q(x) and the net outward flux, W/m2)

Ground truth: the Kirchhoff transform theta = k0 (T + alpha T^2 / 2) turns the problem into theta'' = -q, which has an exact solution. The verifier compares the solver output against it.
