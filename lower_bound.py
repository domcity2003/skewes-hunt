"""Rule out a crossover in a short interval above the published count at 10**20."""

from mpmath import iv


START = 10**20
PRIME_COUNT_AT_START = 2_220_819_602_560_918_840
WIDTH = 1_600_000_000
SIEVE_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)


def count_survivors(start, stop):
    """Count integers in (start, stop] coprime to the sieve primes."""
    terms = [(1, 1)]
    for prime in SIEVE_PRIMES:
        terms += [(divisor * prime, -sign) for divisor, sign in terms]
    return sum(sign * (stop // divisor - start // divisor) for divisor, sign in terms)


def li_lower(x):
    # Ei(t) = Euler's constant + log(t) + sum(t**k / (k * k!)).
    # Here t > 0, so a finite sum gives a lower bound.
    t = iv.log(x)
    term = iv.mpf(1)
    total = iv.euler + iv.log(t)
    for k in range(1, 161):
        term *= t / k
        total += term / k
    return total.a


def main():
    iv.dps = 50
    stop = START + WIDTH
    survivors = count_survivors(START, stop)
    count_upper = PRIME_COUNT_AT_START + survivors
    integral_lower = li_lower(START)
    if not integral_lower > count_upper:
        raise ArithmeticError("The prime-count upper bound exceeds the target")
    print(f"Excluded interval: [{START:,}, {stop:,}]")
    print(f"Sieve survivors in the interval above its left endpoint: {survivors:,}")
    print(f"Lower bound on li(x) - pi(x) throughout: {int(integral_lower) - count_upper:,}")


if __name__ == "__main__":
    main()
