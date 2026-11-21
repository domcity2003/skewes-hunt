"""Re-evaluate Smith's published sum with the finite-height Revers theorem."""

from decimal import Decimal
from fractions import Fraction

from mpmath import iv

from replay_smith import ALPHA, CENTER, CUTOFF, ZERO_SUM_LOWER


VERIFIED_HEIGHT = 200_000_000_000
PUBLISHED_FINITE_HEIGHT = 3_000_000_000_000
HALF_WIDTH = "0.00000000585"


def revers_errors(eta):
    """Enclose R1 through R6 of Revers, arXiv:2501.04488v2, Theorem 3.1."""
    b = iv.mpf
    w, eta, alpha, t, a = map(b, (CENTER, eta, ALPHA, CUTOFF, VERIFIED_HEIGHT))
    low, high = w - eta, w + eta
    if not (low > b("44.22") and 0 < eta < w / 100 and 5 * a / (4 * w) < alpha < a**2
            and 2 * iv.pi * iv.exp(1) < t < a < PUBLISHED_FINITE_HEIGHT):
        raise ValueError("Revers's hypotheses failed")
    root = iv.sqrt(alpha)
    log2 = iv.log(2)
    gaussian = iv.exp(-alpha * eta**2 / 2)
    reciprocal = 1 / eta
    if not reciprocal < b("0.082") * alpha:
        raise ValueError("Unexpected minimum in R2")
    r1 = 2 / low + 8 / low**2 + b("58.56") / low**3
    r1 += log2 * high * iv.exp(-low / 2) + 2 * high * iv.exp(-low / 6) / log2
    r2 = b("0.037") / (w * root) * (reciprocal * gaussian + (1 - gaussian) / low)
    r3 = b("0.074") * root * gaussian
    r4 = iv.exp(-t**2 / (2 * alpha)) * (1 + 1 / (w * t)) * (
        alpha / (iv.pi * t**2) * iv.log(t / (2 * iv.pi)) + 8 * iv.log(t) / t + 4 * alpha / t**3
    )
    r5 = b("0.003") / low**2
    r6 = (1 + 22 / (a * w)) * a * iv.log(a) * (
        b("8.283") / a * iv.exp(-alpha * eta**2 / 4 + high / 2)
        + b("7.152") * eta * iv.exp(-a**2 / (2 * alpha) + high / 2)
    )
    return r1, r2, r3, r4, r5, r6


def rational_budget():
    """An independent, deliberately looser check using exact arithmetic."""
    q = Fraction
    w, eta, alpha, t, a = map(q, (CENTER, HALF_WIDTH, ALPHA, CUTOFF, VERIFIED_HEIGHT))
    low, high = w - eta, w + eta

    def exp_exceeds(argument, bound):
        # All omitted terms of this Taylor series are positive.
        term = total = q(1)
        for k in range(1, 81):
            term *= q(argument) / k
            total += term
        if total <= bound:
            raise ArithmeticError(f"Failed exp({argument}) > {bound}")

    def require(condition):
        if not condition:
            raise ArithmeticError("An exact budget inequality failed")

    require(727 < low < high < 728 and t > 24 and t < a <= PUBLISHED_FINITE_HEIGHT)
    require(5 * a / (4 * w) < alpha < a**2 and 0 < eta < w / 100)
    require((6 * 10**9)**2 < alpha < (7 * 10**9)**2)
    require(alpha * eta**2 / 2 > 700 and 1 / eta < q("0.082") * alpha)
    require(t**2 / (2 * alpha) > q("10.4"))
    require(alpha * eta**2 / 4 - high / 2 > 21)
    require(a**2 / (2 * alpha) - high / 2 > 80)
    exp_exceeds("10.4", 32_000)
    exp_exceeds(23, t / 6)
    exp_exceeds(25, t)
    exp_exceeds(27, a)
    exp_exceeds(21, 10**9)

    # Use 3 < pi < 4, e < 3, and 1/2 < log(2) < 1.
    r1 = 2 / low + 8 / low**2 + q("58.56") / low**3 + 5 * high / 2**120
    r2 = q("0.037") / (w * 6 * 10**9) * (q(1, 2**700) / eta + 1 / low)
    r3 = q("0.074") * 7 * 10**9 / 2**700
    r4 = (1 + 1 / (w * t)) / 32_000 * (23 * alpha / (3 * t**2) + 200 / t + 4 * alpha / t**3)
    r5 = q("0.003") / low**2
    r6 = (1 + 22 / (a * w)) * 27 * (q("8.283") / 10**9 + q("7.152") * a * eta / 2**80)
    ceiling = sum((r1, r2, r3, r4, r5, r6))
    require(ceiling < q("0.002775"))
    require(q(ZERO_SUM_LOWER) - 1 - ceiling > q("0.0001"))
    return ceiling


def main():
    iv.dps = 50
    errors = revers_errors(HALF_WIDTH)
    margin = iv.mpf(ZERO_SUM_LOWER) - 1 - sum(errors)
    if not margin > iv.mpf("0.00009"):
        raise ArithmeticError("The conservative margin did not pass")
    lower = Decimal(CENTER) - Decimal(HALF_WIDTH)
    upper = Decimal(CENTER) + Decimal(HALF_WIDTH)
    print(f"log(x) in [{lower}, {upper}]")
    print(f"cited finite verification height A={VERIFIED_HEIGHT}")
    print(f"cited lower bound for negative zero sum={ZERO_SUM_LOWER}")
    for number, error in enumerate(errors, 1):
        print(f"R{number} = {error}")
    print(f"margin lower bound = {margin}")
    rational_budget()
    print("Exact rational cross-check: sum of errors < 0.002775 and margin > 0.0001")


if __name__ == "__main__":
    main()
