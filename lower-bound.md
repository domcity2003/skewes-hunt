# Ruling out an earlier crossover

**Dominic Bavuso** · November 2026

To identify the first crossover, an existence bound must be paired with a proof that $\pi(x)\leq\mathrm{li}(x)$ everywhere earlier. Büthe proved this inequality through $10^{19}$. Axler cites a result of Dusart extending it through $10^{20}$. The calculation below starts at $10^{20}$ and establishes the inequality on a further interval using a published prime count. [Büthe, 2018](https://doi.org/10.1090/MCOM/3264) · [Axler, Remark 7](https://math.colgate.edu/~integers/y34/y34.pdf) · [Dusart, Lemma 2.2](https://doi.org/10.1007/s11139-017-9984-4)

## A finite extension

Let $a=10^{20}$ and $b=a+1.6\times 10^9$. Let $P$ be the product of the first 15 primes, ending at 47. Every prime in $(a,b]$ is coprime to $P$. Inclusion–exclusion gives the number of integers in that interval that survive the sieve.

$$
S(a,b)
=\sum_{d\mid P}\mu(d)
\left(\left\lfloor\frac bd\right\rfloor-\left\lfloor\frac ad\right\rfloor\right)
=221{,}926{,}542.
$$

Using Deléglise’s tabulated $\pi(a)=2{,}220{,}819{,}602{,}560{,}918{,}840$, monotonicity gives, for every $x\in[a,b]$,

$$
\mathrm{li}(x)-\pi(x)
\geq \mathrm{li}(a)-\pi(a)-S(a,b)
>818{,}101.
$$

The [replay](lower_bound.py) evaluates the sieve exactly and lower-bounds $\mathrm{li}(a)$ by truncating the positive-term series for $\mathrm{Ei}(\log a)$ with interval arithmetic. Run `python lower_bound.py` with mpmath installed. [Deléglise’s prime count](https://t5k.org/notes/md6-96.html) · [NIST DLMF, §6.6](https://dlmf.nist.gov/6.6.E1)

## What remains

The extension is small relative to the gap up to the [existence interval](tapered-bound.md) near $1.397\times10^{316}$. The sieve uses a prime count at one endpoint and stops paying off as the interval grows. A large advance would require a tighter upper bound for $\pi(x)-\mathrm{li}(x)$ over whole intervals, with an error bound valid between sample points. Büthe’s method of bounding prime-counting errors over ranges is a starting point. [Büthe, 2018](https://doi.org/10.1090/MCOM/3264)
