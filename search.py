"""Floating-point experiments for Revers's crossover search (arXiv:2501.04488v2)."""

import argparse
import json
import math
from pathlib import Path

import numpy as np


SMITH_ENDPOINT = 727.9513354023
BENCHMARK_ZERO_COUNT = 2_000_000
BENCHMARK_SUM = -1.006553478788955
BENCHMARK_OMEGA = 727.952018
BENCHMARK_ALPHA = 1.34e11
TABLE5_CENTERS = (41.6522, 84.7316, 136.0262, 154.9746, 157.8305, 175.9619,
                  179.0999, 190.1264, 214.2382, 259.9694, 275.7852, 298.0048,
                  314.0808, 316.1456)


class ZeroSum:
    def __init__(self, ordinates, alpha=None):
        self.ordinates = ordinates
        self.first = 2 / (0.5 + 1j * ordinates)
        self.second = None
        if alpha is not None:
            rho = 0.5 + 1j * ordinates
            self.first *= np.exp(-ordinates**2 / (2 * alpha))
            self.second = self.first / rho

    def evaluate(self, omega, origin=None):
        if origin is None:
            phase = self.ordinates * omega
        else:
            anchor, phase_at_anchor = origin
            phase = self.ordinates * (omega - anchor) + phase_at_anchor
        coefficient = self.first
        if self.second is not None:
            coefficient = coefficient + self.second / omega
        terms = coefficient.real * np.cos(phase) - coefficient.imag * np.sin(phase)
        return float(np.sum(terms))

    def phase_origin(self, omega):
        # Reuse phase reduction for nearby exploratory queries. The benchmark
        # uses the direct expression without an anchor.
        return omega, np.remainder(self.ordinates * omega, math.tau)


def error_terms(omega, eta, alpha, cutoff, verified_height=None):
    """Theorem 3.1. With no verified height supplied, assume full RH."""
    lower, upper = omega - eta, omega + eta
    if eta <= 0 or lower < 44.22 or alpha <= 0 or cutoff <= 2 * math.pi * math.e:
        raise ValueError("Parameters fall outside Theorem 3.1")
    root_alpha = math.sqrt(alpha)
    gaussian = math.exp(-alpha * eta**2 / 2)
    log_two = math.log(2)
    errors = {
        "R1": (2 / lower + 8 / lower**2 + 58.56 / lower**3
               + log_two * upper * math.exp(-lower / 2)
               + 2 * upper * math.exp(-lower / 6) / log_two),
        "R2": (0.037 / (omega * root_alpha)
               * (min(0.082 * alpha, 1 / eta) * gaussian + (1 - gaussian) / lower)),
        "R3": 0.074 * root_alpha * gaussian,
        "R4": (math.exp(-cutoff**2 / (2 * alpha)) * (1 + 1 / (omega * cutoff))
               * (alpha / (math.pi * cutoff**2) * math.log(cutoff / (2 * math.pi))
                  + 8 * math.log(cutoff) / cutoff + 4 * alpha / cutoff**3)),
        "R5": 0.003 / lower**2,
    }
    if verified_height is not None:
        height = verified_height
        if not (cutoff <= height and 5 * height / (4 * omega) <= alpha <= height**2
                and eta < omega / 100):
            raise ValueError("Parameters fall outside the finite-height RH conditions")
        errors["R6"] = ((1 + 22 / (height * omega)) * height * math.log(height)
                        * (8.283 / height * math.exp(-alpha * eta**2 / 4 + upper / 2)
                           + 7.152 * eta * math.exp(-height**2 / (2 * alpha) + upper / 2)))
    return errors


def load_zeros(path):
    ordinates = np.load(path, mmap_mode="r") if path.suffix == ".npy" else np.loadtxt(path)
    if (ordinates.ndim != 1 or len(ordinates) < BENCHMARK_ZERO_COUNT
            or not math.isclose(ordinates[0], 14.134725142, rel_tol=0, abs_tol=5e-9)
            or not np.all(np.isfinite(ordinates)) or not np.all(np.diff(ordinates) > 0)):
        raise ValueError("Expected an ordered table starting at the first zeta zero")
    return ordinates


def benchmark(ordinates):
    cutoff = 1_131_944.4719
    zeros = ordinates[:BENCHMARK_ZERO_COUNT]
    if not zeros[-1] < cutoff < ordinates[BENCHMARK_ZERO_COUNT]:
        raise ValueError("The table does not bracket the published two-million-zero cutoff")
    zero_sum = ZeroSum(zeros, BENCHMARK_ALPHA).evaluate(BENCHMARK_OMEGA)
    # The source ordinates are accurate to 4e-9 and the sum here uses float64.
    if abs(zero_sum - BENCHMARK_SUM) > 5e-7:
        raise ArithmeticError("Failed the Chao–Plymen benchmark")
    margins = {}
    for eta in (0.00016, 0.0001061, 0.0001060):
        errors = error_terms(BENCHMARK_OMEGA, eta, BENCHMARK_ALPHA, cutoff, 1.022e7)
        margins[str(eta)] = -1 - zero_sum - math.fsum(errors.values())
    if not margins["0.0001061"] > 0 > margins["0.000106"]:
        raise ArithmeticError("Failed to reproduce the benchmark margin sign change")
    return {"zero_sum": zero_sum, "difference_from_revers": zero_sum - BENCHMARK_SUM,
            "margins_using_finite_height_rh": margins}


def scan_early_regions(ordinates):
    coarse = ZeroSum(ordinates[:100_000])
    fine = ZeroSum(ordinates[:1_000_000])
    log_ten = math.log(10)
    results = []
    for center in TABLE5_CENTERS:
        coarse_origin = coarse.phase_origin(center * log_ten)
        fine_origin = fine.phase_origin(center * log_ten)
        grid = np.linspace(center - 0.0005, center + 0.0005, 151)
        scores = [-1 - coarse.evaluate(point * log_ten, coarse_origin) for point in grid]
        shortlist = grid[np.argsort(scores)[-8:]]
        # Recheck nearby points at one million zeros, plus the printed center.
        local_grid = np.unique(np.r_[center, (shortlist[:, None] + np.linspace(-8e-6, 8e-6, 21)).ravel()])
        best_score, best_point = max((-1 - fine.evaluate(point * log_ten, fine_origin), point)
                                     for point in local_grid)
        results.append({"table5_log10_x": center, "best_sample_log10_x": best_point,
                        "score_at_one_million": best_score})
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("zeros", type=Path, help="Odlyzko text/gzip or a NumPy ordinate array")
    parser.add_argument("--early", action="store_true", help="scan the neighborhoods in Revers Table 5")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    ordinates = load_zeros(args.zeros)
    report = {"source": args.zeros.name, "available_zeros": len(ordinates),
              "arithmetic": "float64 exploratory", "benchmark": benchmark(ordinates)}
    if args.early:
        report["table5_scan"] = scan_early_regions(ordinates)
    output = json.dumps(report, indent=2) + "\n"
    if args.output is None:
        print(output, end="")
    else:
        args.output.write_text(output, encoding="utf-8")


if __name__ == "__main__":
    main()
