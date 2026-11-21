# The Skewes crossover

**Dominic Bavuso** · November 2026

How early can we prove that the number of primes up to *x* exceeds its standard estimate, li(*x*)?

**Candidate result.** Using Smith and Platt’s published zero sum, the tapered Gaussian argument places a crossover in this range:

> **727.95133539940 ≤ log *x* ≤ 727.95133540208**

The right endpoint is **0.22 billionths earlier on the log scale** than the one in Smith’s 2016 thesis. Each bar shows an interval with at least one crossover. Its right endpoint gives a latest possible location for the first crossover. [Smith, Theorem 9.1](https://etheses.whiterose.ac.uk/id/eprint/16409/)

![On the log scale, Smith’s interval ends at an offset of +1.30 billionths and the tapered Gaussian interval ends at +1.08 billionths.](figures/crossover-intervals.png)

[Derivation](tapered-bound.md) · [Quick replay](#quick-replay) · [Full numerical audit](audit.md) · [Search notes](experiments.md)

## What’s verified

The arithmetic replay passes with the published sum. The independent calculation of that sum from **103.8 billion zeta zeros** is underway. The derivation and literature priority still need external review. The exact first crossover remains open.

The method takes a weighted average of the gap between the prime count and its estimate. The weight fades smoothly to zero at the ends of the interval, which lets the argument bound the effect of the highest zeros. A positive weighted average proves the gap is positive somewhere in the interval. [Derivation and taper figure](tapered-bound.md) · [Platt and Trudgian’s finite verification](https://arxiv.org/abs/2004.09765)

## Quick replay

Run with Python 3.12. This checks the error budget against the published zero sum without downloading the archive. The positive margin exceeds **0.0000326**.

~~~sh
python -m pip install mpmath
python tapered.py
~~~

| Further reading | What it covers |
| :--- | :--- |
| [Independent numerical audit](audit.md) | Input checks, rounding bounds, and a resumable replay of the full archive |
| [Search experiments](experiments.md) | Tests of an earlier candidate using 158.2 million zeros |
| [Lower-bound calculation](lower-bound.md) | A sieve that rules out a short interval just above 10<sup>20</sup> |
