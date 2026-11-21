"""Sweep Gaussian widths at the earlier Saouter–Trudgian–Demichel candidate."""

import argparse
import json
from pathlib import Path

import numpy as np

from replay_smith import C_LOWER, smith_errors


TARGET = 727.951332982
CONTROL = 727.951335792
ENDPOINT = 727.95133540175
ALPHAS = np.sort(np.r_[np.geomspace(5e11, 4e13, 25), 6e12])
CUTOFFS = (1_131_944.4719, 2_546_000, 4_646_000, 6_746_000, 8_846_000, 10_379_600)


def accumulate(total, ordinates, centers, alphas):
    for begin in range(0, len(ordinates), 250_000):
        batch = ordinates[begin:begin + 250_000]
        rho = 0.5 + 1j * batch
        phase = np.multiply.outer(centers - TARGET, batch)
        phase += np.remainder(batch * TARGET, 2 * np.pi)
        coefficient = 2 / rho + np.multiply.outer(1 / centers, 2 / rho**2)
        terms = coefficient.real * np.cos(phase) - coefficient.imag * np.sin(phase)
        weights = np.exp(np.multiply.outer(batch**2, -0.5 / alphas))
        total -= terms @ weights


def scan(zeros):
    centers = np.r_[TARGET + np.linspace(-3e-7, 3e-7, 25), CONTROL]
    budgets = np.maximum(ENDPOINT - centers - 1e-12, 4e-6 * (centers == CONTROL))
    total = np.zeros((len(centers), len(ALPHAS)))
    results = []
    start = 0
    for cutoff in CUTOFFS:
        if cutoff >= zeros[-1]:
            break
        stop = np.searchsorted(zeros, cutoff, side="right")
        accumulate(total, zeros[start:stop], centers, ALPHAS)
        start = stop
        if cutoff == 10_379_600:
            control = total[-1, np.searchsorted(ALPHAS, 6e12)]
            if stop != 22_000_000 or abs(control - 1.002889274) > 5e-8:
                raise ArithmeticError(f"Failed the Teal Sea 22-million-zero control: {control}")
        rows = []
        for ci, (center, eta) in enumerate(zip(centers, budgets)):
            for ai, alpha in enumerate(ALPHAS):
                errors = smith_errors(str(eta), str(center), str(alpha), str(cutoff))
                error = sum(float(term.b) for term in errors)
                score = total[ci, ai]
                rows.append({"center": center, "eta": eta, "alpha": alpha, "sum": score,
                             "error": error, "score_minus_main_term": score + float(C_LOWER),
                             "margin_estimate": score + float(C_LOWER) - error})
        earlier = [row for row in rows if row["center"] != CONTROL]
        rank = lambda row: row["margin_estimate"]
        results.append({"cutoff": cutoff, "zeros": int(stop),
                        "best_at_earlier_target": max(earlier, key=rank),
                        "highest_raw_signal_at_target": max(earlier, key=lambda row: row["sum"]),
                        "best_at_control": max((row for row in rows if row["center"] == CONTROL), key=rank)})
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("zeros", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    zeros = np.load(args.zeros, mmap_mode="r", allow_pickle=False)
    report = {"arithmetic": "float64 zero sums, screening only, assuming RH",
              "target_center": TARGET, "tested_centers": 25, "tested_alphas": ALPHAS.tolist(),
              "results": scan(zeros)}
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
