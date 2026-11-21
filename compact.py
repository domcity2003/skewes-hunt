"""Probe a compactly supported alternative to the Gaussian Lehman weight."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import spherical_jn

from earlier import TARGET, ENDPOINT


CUTOFF = 10_379_600
ETA = 2.1e-6
ORDERS = tuple(range(3, 11))


def double_factorial(p):
    return math.prod(range(1, 2 * p + 2, 2))


def zero_tail(order, cutoff):
    """Saouter–Trudgian–Demichel Lemma 2.7 with phi(t) = t**(-order).

    Returns a bound scaled by cutoff**(order - 1).
    """
    decay = order - 1
    return ((math.log(cutoff / (2 * math.pi)) / decay + 1 / decay**2) / (2 * math.pi)
            + (4 * math.log(cutoff) + 2 / order) / cutoff)


def compact_tail(p, eta, cutoff):
    # Apply the finite spherical Hankel expansion to bound |j_p| termwise.
    return 2 * double_factorial(p) * sum(
        math.factorial(p + k) / (2**k * math.factorial(k) * math.factorial(p - k))
        * zero_tail(p + k + 2, cutoff) / (eta * cutoff)**(p + k + 1)
        for k in range(p + 1)
    )


def scan(zeros, cutoff):
    centers = TARGET + np.linspace(-3e-7, 3e-7, 13)
    if centers[-1] + ETA >= ENDPOINT:
        raise ValueError("A window exceeds the endpoint")
    count = int(np.searchsorted(zeros, cutoff, side="right"))
    totals = np.zeros((len(centers), len(ORDERS)))
    for start in range(0, count, 250_000):
        gamma = zeros[start:min(start + 250_000, count)]
        reduced = np.remainder(gamma * TARGET, math.tau)
        phase = reduced + np.multiply.outer(centers - TARGET, gamma)
        rho = 0.5 + 1j * gamma
        coefficient = 2 / rho + np.multiply.outer(1 / centers, 2 / rho**2)
        head = coefficient.real * np.cos(phase) - coefficient.imag * np.sin(phase)
        argument = ETA * gamma
        windows = np.column_stack([
            double_factorial(p) * spherical_jn(p, argument) / argument**p for p in ORDERS
        ])
        totals -= head @ windows
    report = []
    for index, p in enumerate(ORDERS):
        best = int(np.argmax(totals[:, index]))
        center, score = centers[best], totals[best, index]
        low, high = center - ETA, center + ETA
        prime_cost = (1 + 2 / low + 8 / low**2 + 58.56 / low**3
                      + high * math.log(2) * math.exp(-low / 2)
                      + 2 * high / math.log(2) * math.exp(-low / 6))
        expansion_cost = 0.00292 / low**2 + 0.046211 * ETA / (low * center)
        remainder = compact_tail(p, ETA, cutoff) * (1 + 1 / (center * cutoff))
        report.append({"order": p, "center": center, "zero_sum": score,
                       "tail_allowance": remainder, "prime_cost": prime_cost,
                       "margin_estimate": score - prime_cost - expansion_cost - remainder})
    return {"cutoff": cutoff, "zero_count": count, "weight_half_width": ETA,
            "center_count": len(centers), "assumption": "RH", "arithmetic": "float64 screening",
            "results": report}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("zeros", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    zeros = np.load(args.zeros, mmap_mode="r", allow_pickle=False)
    if zeros[-1] <= CUTOFF:
        raise ValueError("The data do not cover the requested cutoff")
    report = scan(zeros, CUTOFF)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
