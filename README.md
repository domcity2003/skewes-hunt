# An upper bound for the first prime-counting crossover

**Dominic Bavuso** · November 2026

Let $\pi(x)$ count the primes up to $x$, and let $\operatorname{li}(x)$ be the logarithmic integral, an estimate for that count. The first crossover is the smallest $x$ for which $\pi(x)>\operatorname{li}(x)$.

## Main result

Using the published zero sum of Smith and Platt, the calculation establishes that $\pi(x)>\operatorname{li}(x)$ for some $x$ satisfying

$$
\boxed{727.95133539940\;\leq\;\log x\;\leq\;727.95133540208.}
$$

The upper endpoint is $2.2\times10^{-10}$ lower on the log scale than Smith’s thesis endpoint. The argument uses Platt and Trudgian’s verification of zeta zeros through height $3\times10^{12}$ and bounds the contribution from higher zeros directly. [Smith, Theorem 9.1](https://etheses.whiterose.ac.uk/id/eprint/16409/) · [Platt and Trudgian, Theorem 1](https://arxiv.org/abs/2004.09765)

![Two existence intervals on the log scale. The tapered Gaussian interval ends 0.22 billionths earlier than Smith’s.](figures/crossover-intervals.png)

*Each bar marks an interval containing a point where the prime count exceeds the estimate. The right endpoint therefore gives an upper bound for the first such point.*

## Method and verification

The calculation averages the difference between the actual and estimated prime counts after scaling it by $\log(x)/\sqrt{x}$. A nonnegative weight tapers to zero at the two ends of the interval. This preserves most of the published Gaussian sum and limits the contribution of zeros above the verified height. The resulting weighted integral has a certified lower bound of $3.26\times10^{-5}$ using the cited numerical input. [Research note and derivation](tapered-bound.md) · [Revers’s explicit formula](https://arxiv.org/html/2501.04488v2#S4)

Run the interval check with Python 3.12 and `mpmath`. The published zero sum is supplied as a constant, so the replay requires no zero-table downloads.

```sh
python -m pip install mpmath
python tapered.py
```

The arithmetic has been checked with interval arithmetic and exact rational bounds. Independent review of the derivation, the original 103.8-billion-zero computation, and literature priority remains open. The exact first crossover is unknown.

The [independent replay](audit.md) checks the original zero sum from the LMFDB files and can resume at any completed file.

## Finding the first crossover

An upper bound gives a deadline by which a crossover has occurred. Identifying the first crossover also requires ruling out every earlier value. A finite sieve check excludes a short interval immediately above $10^{20}$, the lower endpoint reached by existing work. The gap between this range and the upper bound near $10^{316}$ remains open. [Lower-bound calculation](lower-bound.md)

## Next experiment

The earlier candidate at $\log x\approx727.951332982$ remains unresolved. A search over 158.2 million zeros recovered a known positive control, while the earlier signal stayed below the threshold in the tested configurations. Smith reports a stronger signal from roughly 100 billion zeros, and its raw sum and parameters are still needed for a replay. [Search results and reproduction commands](experiments.md) · [Smith, Section 10](https://etheses.whiterose.ac.uk/16409/1/The%20hunt%20for%20Skewes%27%20number.pdf)
