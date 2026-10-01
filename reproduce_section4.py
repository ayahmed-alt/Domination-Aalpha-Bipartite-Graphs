"""
Reproducibility code for:
"Domination Number and A_alpha-Spectral Properties of Bipartite Graphs"

This script reproduces the computational consistency checks described in
Section 4 of the manuscript:
  1. exhaustive tests on connected bipartite graphs in the NetworkX Graph Atlas
     up to order 7;
  2. the reproducible random balanced-bipartite experiment with seed 20260926;
  3. selected comparisons at alpha = 0.3;
  4. the n = 20 complete-bipartite phase-transition values.

Requirements:
    Python 3.10+
    numpy
    networkx

The computations are supplementary checks only and are not used in the proofs.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
import numpy as np
import networkx as nx


ALPHAS = (0.0, 0.1, 0.3, 0.5, 0.51, 0.7, 0.9, 0.99)
SEED = 20260926
TOL = 1e-9


def domination_number(G: nx.Graph) -> int:
    """Compute gamma(G) exactly by exhaustive subset enumeration."""
    nodes = list(G.nodes())
    n = len(nodes)
    closed = {v: set(G.neighbors(v)) | {v} for v in nodes}

    for k in range(1, n + 1):
        for S in itertools.combinations(nodes, k):
            dominated = set()
            for v in S:
                dominated.update(closed[v])
            if len(dominated) == n:
                return k
    return n


def aalpha_matrix(G: nx.Graph, alpha: float) -> np.ndarray:
    """Return A_alpha(G) = alpha D(G) + (1-alpha) A(G)."""
    nodes = list(G.nodes())
    A = nx.to_numpy_array(G, nodelist=nodes, dtype=float)
    degrees = np.array([G.degree(v) for v in nodes], dtype=float)
    return alpha * np.diag(degrees) + (1.0 - alpha) * A


def rho_alpha(G: nx.Graph, alpha: float) -> float:
    """Largest eigenvalue of the real symmetric A_alpha matrix."""
    vals = np.linalg.eigvalsh(aalpha_matrix(G, alpha))
    return float(vals[-1])


def theorem31_bound(G: nx.Graph, alpha: float, gamma: int | None = None) -> float:
    """
    Theorem 3.1 bound.

    S = n - gamma + 2.
    For alpha <= 1/2: S/2.
    For alpha > 1/2:
      [alpha S + sqrt(alpha^2 S^2
       - 4(2 alpha - 1) delta (S-delta))]/2.
    """
    n = G.number_of_nodes()
    if gamma is None:
        gamma = domination_number(G)
    delta = min(dict(G.degree()).values())
    S = n - gamma + 2

    if alpha <= 0.5:
        return S / 2.0

    rad = (
        alpha**2 * S**2
        - 4.0 * (2.0 * alpha - 1.0) * delta * (S - delta)
    )
    # Protect only against tiny negative floating-point roundoff.
    if rad < 0 and rad > -1e-12:
        rad = 0.0
    return (alpha * S + np.sqrt(rad)) / 2.0


def complete_bipartite_rho(p: int, q: int, alpha: float) -> float:
    """Closed form for rho_alpha(K_{p,q})."""
    n = p + q
    rad = alpha**2 * n**2 - 4.0 * (2.0 * alpha - 1.0) * p * q
    return (alpha * n + np.sqrt(rad)) / 2.0


def exhaustive_atlas_experiment():
    """Reproduce the exhaustive Graph Atlas counts and Theorem 3.1 checks."""
    counts = {}
    gamma2_counts = {}
    tests = 0
    min_gap = float("inf")
    violations = []

    for G0 in nx.graph_atlas_g():
        n = G0.number_of_nodes()
        if n < 2 or n > 7:
            continue
        if not nx.is_connected(G0) or not nx.is_bipartite(G0):
            continue

        G = nx.convert_node_labels_to_integers(G0)
        gamma = domination_number(G)
        counts[n] = counts.get(n, 0) + 1
        if gamma == 2:
            gamma2_counts[n] = gamma2_counts.get(n, 0) + 1

        for alpha in ALPHAS:
            rho = rho_alpha(G, alpha)
            bound = theorem31_bound(G, alpha, gamma)
            gap = bound - rho
            min_gap = min(min_gap, gap)
            tests += 1
            if gap < -TOL:
                violations.append((n, gamma, alpha, rho, bound, gap))

    return counts, gamma2_counts, tests, min_gap, violations


def random_balanced_bipartite_graph(
    n: int, prob: float, rng: np.random.Generator
) -> nx.Graph:
    """Generate one balanced bipartite Bernoulli graph using the supplied RNG."""
    p = n // 2
    q = n - p
    G = nx.Graph()
    X = list(range(p))
    Y = list(range(p, p + q))
    G.add_nodes_from(X, bipartite=0)
    G.add_nodes_from(Y, bipartite=1)

    draws = rng.random((p, q))
    for i in range(p):
        for j in range(q):
            if draws[i, j] < prob:
                G.add_edge(X[i], Y[j])
    return G


def random_experiment():
    """
    Reproduce the random experiment exactly from the manuscript protocol:
      seed = 20260926
      n in {10,12,14,16}
      edge probabilities in {0.15,0.30,0.50,0.70}
      40 generated graphs per (n, probability)
      disconnected graphs discarded
      eight alpha values tested.
    """
    rng = np.random.default_rng(SEED)
    ns = (10, 12, 14, 16)
    probs = (0.15, 0.30, 0.50, 0.70)

    retained_by_n = {n: 0 for n in ns}
    min_gap_by_n = {n: float("inf") for n in ns}
    violations = []
    tests = 0

    for n in ns:
        for prob in probs:
            for _ in range(40):
                G = random_balanced_bipartite_graph(n, prob, rng)
                if not nx.is_connected(G):
                    continue

                retained_by_n[n] += 1
                gamma = domination_number(G)

                for alpha in ALPHAS:
                    rho = rho_alpha(G, alpha)
                    bound = theorem31_bound(G, alpha, gamma)
                    gap = bound - rho
                    min_gap_by_n[n] = min(min_gap_by_n[n], gap)
                    tests += 1
                    if gap < -TOL:
                        violations.append(
                            (n, prob, gamma, alpha, rho, bound, gap)
                        )

    return retained_by_n, min_gap_by_n, tests, violations


def selected_comparison():
    """Values reported in the manuscript at alpha = 0.3."""
    alpha = 0.3
    graphs = [
        ("P6", nx.path_graph(6)),
        ("C8", nx.cycle_graph(8)),
        ("K3,3", nx.complete_bipartite_graph(3, 3)),
    ]

    rows = []
    for name, G in graphs:
        gamma = domination_number(G)
        rho = rho_alpha(G, alpha)
        bip_bound = theorem31_bound(G, alpha, gamma)
        general_bound = G.number_of_nodes() - gamma
        rows.append((name, G.number_of_nodes(), gamma, rho, bip_bound, general_bound))
    return rows


def phase_transition_table():
    """Values for K_{2,18}, K_{3,17}, K_{5,15}, K_{10,10} at n=20."""
    pairs = ((2, 18), (3, 17), (5, 15), (10, 10))
    alphas = (0.00, 0.25, 0.49, 0.50, 0.51, 0.75, 0.90)
    rows = []
    for alpha in alphas:
        rows.append(
            (alpha,) + tuple(complete_bipartite_rho(p, q, alpha) for p, q in pairs)
        )
    return rows


def main():
    print("=== Exhaustive Graph Atlas experiment ===")
    counts, gamma2, tests, min_gap, violations = exhaustive_atlas_experiment()
    print("Connected bipartite graph counts:", counts)
    print("gamma=2 counts:", gamma2)
    print("Spectral tests:", tests)
    print(f"Minimum bound-minus-radius gap: {min_gap:.16g}")
    print("Violations:", len(violations))

    print("\n=== Random experiment ===")
    retained, min_gaps, random_tests, random_violations = random_experiment()
    print("Seed:", SEED)
    print("Retained connected graphs:", retained)
    print("Spectral tests:", random_tests)
    print("Minimum gaps by n:")
    for n in sorted(min_gaps):
        print(f"  n={n}: {min_gaps[n]:.16g}")
    print("Violations:", len(random_violations))

    print("\nCombined spectral tests:", tests + random_tests)

    print("\n=== Selected comparison at alpha=0.3 ===")
    print("Graph   n  gamma   rho_alpha   bipartite_bound   general_bound")
    for row in selected_comparison():
        name, n, gamma, rho, bb, gb = row
        print(f"{name:5s} {n:2d} {gamma:5d} {rho:11.6f} {bb:17.6f} {gb:14.6f}")

    print("\n=== Phase transition, n=20 ===")
    print("alpha   K2,18       K3,17       K5,15       K10,10")
    for row in phase_transition_table():
        alpha, *vals = row
        print(f"{alpha:4.2f}  " + "  ".join(f"{v:10.6f}" for v in vals))


if __name__ == "__main__":
    main()
