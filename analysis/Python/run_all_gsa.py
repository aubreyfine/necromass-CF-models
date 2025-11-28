"""
run_all_gsa.py

Convenience driver script to run all global sensitivity analyses and
rB significance diagnostics for the necromass-CF framework.

This script:
  1. Runs Sobol, Delta, and Morris analyses (gsa_numeric.py)
  2. Runs rB significance profiling (rb_significance.py)

Outputs are written to:
  data/derived/gsa/

Author: A.K. Fine
"""

import os

# Ensure we are in the repo root (so relative paths work)
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(REPO_ROOT)

from analysis.Python import gsa_numeric, rb_significance  # type: ignore


def main():
    print("=== necromass-CF: global sensitivity workflow ===")
    print("Working directory:", os.getcwd())
    print("\n[1/2] Running Sobol, Delta, and Morris analyses...")
    gsa_numeric.main(N=4096)

    print("\n[2/2] Running rB significance analysis...")
    rb_significance.rB_significance_analysis()

    print("\nAll analyses complete.")
    print("Tidy numeric outputs written to: data/derived/gsa/")
    print("These can be visualized using the R scripts in analysis/R/.")


if __name__ == "__main__":
    main()
