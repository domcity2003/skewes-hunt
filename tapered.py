"""Check the finite-height error budget for the tapered Gaussian."""

from decimal import Decimal
import json
from pathlib import Path

from mpmath import iv

from replay_smith import ALPHA, CENTER, CUTOFF
from unconditional import PUBLISHED_FINITE_HEIGHT, common_errors


LEFT_WIDTH = "1.60e-9"
RIGHT_WIDTH = "1.08e-9"
TAPER_WIDTH = "0.71e-9"
CAUCHY_RADIUS = "0.42e-9"
DERIVATIVE_ORDER = 390


def error_budget(cutoff=CUTOFF):
    b = iv.mpf
    w, alpha, t, a = map(b, (CENTER, ALPHA, cutoff, PUBLISHED_FINITE_HEIGHT))
    left, right, width, radius = map(b, (LEFT_WIDTH, RIGHT_WIDTH, TAPER_WIDTH, CAUCHY_RADIUS))
    order = DERIVATIVE_ORDER
    boxes = order + 1
    low, high = w - left, w + right
    amplitude = iv.sqrt(alpha / (2 * iv.pi))
    rate = 2 * boxes / width
    decay = alpha * (1 - alpha * width**2 / (12 * boxes))
    if not (low > b("44.22") and 2 * iv.pi * iv.exp(1) < t < a
            and order >= 1 and 0 < width < min(left, right)
            and 0 < radius < low / 5 and 0 < decay and order < rate * radius):
        raise ValueError("The taper bounds' hypotheses failed")

    mass_loss = sum(amplitude * iv.exp(-decay * h**2 / 2) / (decay * h)
                    for h in (left - width / 2, right - width / 2))
    reciprocal_squares = 2 * b("0.0231050")
    reciprocal_sum = iv.log(a / (2 * iv.pi))**2 / (2 * iv.pi) + 2 * b("0.9321")
    taper_error = mass_loss * (reciprocal_sum + reciprocal_squares / w)
    reciprocal_error = reciprocal_squares * (left + right) / (low * w)
    high_zero_error = (
        2 * iv.log(a) * (left + right) * amplitude * (1 + 4 / (a * (low - radius)))
        * (rate / a)**order * iv.exp(high / 2 + alpha * radius**2 / 2)
        / (1 - order / (rate * radius))
    )
    r1, r4, r5 = common_errors(w, low, high, alpha, t)
    return {
        "prime-side R1": r1,
        "taper mass error": taper_error,
        "freezing 1/u": reciprocal_error,
        "Gaussian zero tail R4": r4,
        "exponential-integral remainder R5": r5,
        "unverified zeros above A": high_zero_error,
    }


def main():
    iv.dps = 60
    result_path = Path(__file__).with_name("results") / "zero-sum.json"
    result = json.loads(result_path.read_text(encoding="utf-8"))
    if result["settings"] != {"center": CENTER, "alpha": ALPHA, "cutoff": CUTOFF}:
        raise ValueError("The released zero sum uses different Gaussian parameters")
    zero_sum_lower = result["enclosure"][0]
    margin_floor = result["margin_lower_bound"]
    errors = error_budget()
    margin = iv.mpf(zero_sum_lower) - 1 - sum(errors.values())
    if not margin > iv.mpf(margin_floor) > 0:
        raise ArithmeticError("The conservative margin did not pass")
    lo = Decimal(CENTER) - Decimal(LEFT_WIDTH)
    hi = Decimal(CENTER) + Decimal(RIGHT_WIDTH)
    print(f"log(x) in [{lo}, {hi}]")
    print(f"cited finite verification height A={PUBLISHED_FINITE_HEIGHT}")
    print(f"audit lower bound for negative zero sum={zero_sum_lower}")
    for label, error in errors.items():
        print(f"{label} majorant: {error}")
    print(f"lower bound from these majorants: {margin}")
    print(f"Conservative margin lower bound: {margin_floor}")


if __name__ == "__main__":
    main()
