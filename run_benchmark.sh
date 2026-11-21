#!/bin/bash
set -euo pipefail

echo "=== Running reference solver ==="
python3 /app/solution/reference_solver.py

echo "=== Running pytest verifier ==="
pytest /app/tests/test_benchmark.py -v
