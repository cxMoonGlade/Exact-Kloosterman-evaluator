#!/usr/bin/env python3
"""Run the unchanged historical Direct PARI worker on a recorded-grid input.

This input-checking launcher was added for publication and was not timed in
2026-09-30's saved campaign. It never writes the archived measurements.
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("N", type=int, choices=(2, 4, 6, 8, 12))
    parser.add_argument("L", type=int, choices=(2, 3))
    parser.add_argument("a", type=int, choices=(1, 2, 3))
    args = parser.parse_args()
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    subprocess.run([sys.executable, str(Path(__file__).with_name("direct_worker.py")),
                    str(args.N), str(args.L), str(args.a)],
                   env=env, check=True, timeout=120)


if __name__ == "__main__":
    main()
