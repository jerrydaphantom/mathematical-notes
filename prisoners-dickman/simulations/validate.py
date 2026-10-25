"""Reproduce the 100-prisoners note: exact probabilities, Monte Carlo, and figures.

Usage from `prisoners-dickman/`:
    python -m pip install -r requirements.txt
    python simulations/validate.py --trials 100000 --seed 20261025

No network access, user credentials, or GitHub integration is needed.
The simulation samples random permutations independently of the exact recurrence.
"""

from __future__ import annotations

import argparse
import csv
from itertools import permutations
from math import sqrt
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch

HERE = Path(__file__).resolve().parents[1]
FIGURES = HERE / 'figures'
CASES = [(100, 50), (100, 40), (100, 30), (150, 50)]


def exact_success(n: int, m: int) -> float:
    """Exact recurrence q(n,m) = (1/n)*sum_{k=1}^{min(n,m)}q(n-k,m).

    Direct summation avoids floating-point cancellation of sliding windows.
    """
    if n < 0 or m < 0:
        raise ValueError('n and m must be nonnegative')
    if n == 0:
        return 1.0
    if m == 0:
        return 0.0
    q = [1.0] + [0.0] * n
    for j in range(1, n+1):
        q[j] = sum(q[j-k] for k in range(1, min(j, m)+1)) / j
    return q[n]


def longest_cycle(perm) -> int:
    """Lengths of cycles of a permutation represented as a 0-based mapping."""
    n = len(perm)
    seen = [False] * n
    longest = 0
    for i in range(n):
        if seen[i]:
            continue
        j, length = i, 0
        while not seen[j]:
            seen[j] = True
            length += 1
            j = int(perm[j])
        longest = max(longest, length)
    return longest


def monte_carlo(n: int, m: int, trials: int, rng: np.random.Generator):
    successes = sum(longest_cycle(rng.permutation(n)) <= m for _ in range(trials))
    observed = successes / trials
    stderr = sqrt(observed*(1-observed) / trials)
    return observed, stderr


def check_small_n(nmax=7):
    """Independent brute force check for all small permutation sizes and limits."""
    for n in range(1, nmax+1):
        lengths = [longest_cycle(p) for p in permutations(range(n))]
        for m in range(0, n+1):
            brute = sum(length <= m for length in lengths) / len(lengths)
            calc = exact_success(n, m)
            assert abs(brute-calc) < 2e-13, (n,m,brute,calc)
    print(f'Brute force cross-check: n=1,...,{nmax}, every m — passed.')


def dickman_curve(u_max=4.3, step=0.0005):
    """Solve rho'(u)=-rho(u-1)/u by trapezoidal method of steps."""
    delay = round(1 / step)
    u = np.arange(round(u_max / step) + 1, dtype=float) * step
    rho = np.ones(len(u), dtype=float)
    for i in range(delay+1, len(u)):
        f0 = rho[i-1-delay] / u[i-1]
        f1 = rho[i-delay] / u[i]
        rho[i] = rho[i-1] - step*(f0+f1)/2
    return u, rho


def decorate(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(alpha=.16, linewidth=.6)
