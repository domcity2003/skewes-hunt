# When the prime count overtakes its estimate

**Dominic Bavuso** · Computational number theory · Research draft

| [Read the paper (PDF)](paper/crossover.pdf) | [Check the calculation](#quick-replay) | [Inspect the results](audit.md#full-replay) |
| :--- | :--- | :--- |

The logarithmic integral li(x) estimates how many primes there are up to x. This project asks how early the actual count, π(x), can be shown to exceed that estimate.

The tapered Gaussian argument gives a point with π(x) > li(x) in the interval

$$
727.95133539940 \le \log x \le 727.95133540208.
$$

The right endpoint is **2.2 × 10⁻¹⁰ earlier on the log scale** than the endpoint in Smith’s 2016 thesis. It gives an upper bound for the first crossover. The argument uses the published verification of zeta zeros through height 3 × 10¹². [Smith, Theorem 9.1](https://etheses.whiterose.ac.uk/id/eprint/16409/) · [Platt and Trudgian, Theorem 1](https://arxiv.org/abs/2004.09765)

![Smith’s interval ends at +1.30 billionths from the center. The tapered Gaussian interval ends at +1.08 billionths.](figures/crossover-intervals.png)

## The calculation

The proof averages the prime-counting error over the interval. Smoothing the weight at both ends gives enough decay to bound the effect of zeros above the verified height. A positive average implies the prime count exceeds its estimate somewhere inside the interval.

| Evidence | Result |
| :--- | :--- |
| [Full replay](audit.md#full-replay) | **103,800,788,359** ordinates across 14,580 files |
| [Two-term sum](results/zero-sum.json) | Enclosure **[1.0028628275787, 1.0028644487398]** |
| [Error budget](tapered-bound.md#arithmetic-check) | Positive margin **> 0.00001495** after all allowances |
| [Lean check](formal/README.md) | Repeated integration by parts and the decay bound for one zero, given the analytic hypotheses |

The replay gives a smaller sum than Smith’s Table 4. The bound above uses the smaller enclosure. External review of the argument and numerical error budget is pending. The exact first crossover is open.

## Quick replay

From a checkout with Python 3.12, run the interval check against the published result file. It needs only mpmath and runs without the zero archive.

```sh
python -m pip install mpmath==1.3.0
python tapered.py
```

Expected output includes `log(x) in [727.95133539940, 727.95133540208]` and `Conservative margin lower bound: 0.00001495`.

The [per-file sums](results/README.md) let readers reconstruct the reported total from the released checkpoint. The [full audit instructions](audit.md#running-the-replay) reproduce the summation from the original LMFDB files and transfer about 1.35 TB. The audit relies on LMFDB for isolation and completeness of the zeros.

## Research files

| File | Contents |
| :--- | :--- |
| [LaTeX source](paper/crossover.tex) | Manuscript and bibliography |
| [Proof notes](tapered-bound.md) | Taper construction, high-zero estimate, and arithmetic |
| [Search notes](experiments.md) | Tests of an earlier candidate |
| [Lower-bound computation](lower-bound.md) | A short interval excluded just above 10²⁰ |

**Suggested citation.** Dominic Bavuso, *A tapered Gaussian bound for the prime-counting crossover*, research draft, November 2026.
