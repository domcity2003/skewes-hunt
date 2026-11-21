"""Test narrower windows using Smith's published Appendix C calculation."""

from decimal import Decimal

from mpmath import iv


CENTER = "727.951335401"
SOURCE_ETA = "1.3e-9"
CUTOFF = 30_610_046_000
ALPHA = "4.5e19"
ZERO_SUM_LOWER = "1.0028804978364"
C_LOWER = "-1.002763"


def smith_errors(eta, center=CENTER, smoothing=ALPHA, height=CUTOFF):
    """Evaluate Smith Theorem 4.1 with M = T."""
    w, eta, alpha, cutoff = map(iv.mpf, (center, eta, smoothing, height))
    lower = w - eta
    scale = alpha * eta
    kernel = iv.sqrt(alpha / (2 * iv.pi)) * iv.exp(-alpha * eta**2 / 2)
    cutoff_log = iv.log(cutoff / (2 * iv.pi))
    cutoff_tail = iv.exp(-cutoff**2 / (2 * alpha))
    quarter_tail = iv.exp(1 / (32 * alpha) - w / 4)
    return (
        iv.mpf("1.812") * (w + eta) * iv.exp(-lower / 6),
        iv.mpf("0.024") * iv.exp(-lower / 4) * (1 + 4 / lower)
        + quarter_tail * (iv.mpf("1.301") + iv.mpf("0.04") * alpha)
        + 2 * quarter_tail / (scale - iv.mpf("0.25"))
        * (iv.log((4 * scale - 1) / (2 * iv.pi))**2 + iv.log(4 * scale - 1) + iv.mpf("0.9321")),
        alpha * cutoff_tail / (2 * iv.pi) * cutoff_log * (1 / cutoff**3 + 2 / cutoff**2),
        kernel / scale * (iv.log(scale / iv.pi)**2 / iv.pi + 4 * iv.log(2 * scale) + iv.mpf("4.52")),
        2 * alpha * cutoff_tail * cutoff_log / (2 * iv.pi * w * cutoff**3)
        + kernel / (w * scale) * (iv.mpf("0.047") + 1 / scale)
        + iv.mpf("0.019") / (iv.sqrt(alpha) * w**2),
        iv.mpf("0.00292") / lower**2,
    )


def main():
    iv.dps = 45
    for eta in (SOURCE_ETA, "8e-10", "7.5e-10", "7e-10"):
        # Smith's C is the negative integral of a positive function here. For
        # smaller eta, the lower bound from the original interval still applies.
        if not 0 < Decimal(eta) <= Decimal(SOURCE_ETA):
            raise ValueError("Window exceeds the source interval")
        error = sum(smith_errors(eta))
        margin = iv.mpf(ZERO_SUM_LOWER) + iv.mpf(C_LOWER) - error
        upper = Decimal(CENTER) + Decimal(eta)
        print(f"eta={eta}  upper log(x)={upper}  error={error}  margin bound={margin}")
    if not sum(smith_errors("7.5e-10")) < iv.mpf("0.000080"):
        raise ArithmeticError("The conservative bound at eta=7.5e-10 did not pass")


if __name__ == "__main__":
    main()
