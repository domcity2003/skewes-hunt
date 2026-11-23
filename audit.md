# Independent replay of Smith and Platt’s zero sum

**Dominic Bavuso** · November 2026

The [tapered Gaussian bound](tapered-bound.md) uses a two-term zero sum recomputed from LMFDB's published ordinates. The [replay](audit_sum.py) checks each file against the archive's checksum and validates block continuity and cumulative counts. The archive supplies isolation and completeness of the zeros. [Smith, Equations (49)–(50) and Table 4](https://etheses.whiterose.ac.uk/id/eprint/16409/) · [LMFDB dataset](https://www.lmfdb.org/zeros/zeta/)

## Full replay

The run finished on November 23, 2026 after reading 14,580 files and 103,800,788,359 ordinates through height 30,610,046,000. It checked about 1.35 TB of source files. The binary64 sum was 1.0028636381592315. The arithmetic allowance gives the outward-rounded enclosure **[1.0028628275787, 1.0028644487398]**. Its lower end leaves more than 0.00001495 in the weighted inequality. [Machine-readable result](results/zero-sum.json) · [Per-file sums](results/file-sums.csv.gz)

Smith's Table 4 value of 1.00288049783641096204 lies outside the enclosure. A diagnostic through height 10,946,000 put the second term's contribution near −0.00001685969, close to the observed difference. The original computation has yet to be inspected.

## Arithmetic

Within each archive block, the ordinates have the form $\gamma=g_0+\delta$, with $g_0$ a binary64 endpoint and $\delta$ a cumulative integer offset in units of $2^{-101}$. The program preserves that integer offset in two limbs until it forms the local phase. It uses the exact decimal $\omega$, the exact dyadic $g_0$, and a fixed decimal approximation to $2\pi$ to reduce $\omega g_0$ by integer arithmetic. It then adds $\omega\delta$ in binary64 and reduces again. This avoids forming a binary64 phase of order $10^{13}$.

The replay uses rounded-to-nearest binary64 arithmetic. The following estimates apply over the enforced block limits, including $\delta<5000$, $\omega<728$, at most 50,000 terms per block, and $|\gamma|\le 30{,}610{,}046{,}000$.

- The fixed decimal period differs from $2\pi$ by less than $10^{-48}$. It is multiplied by fewer than $4\times10^{12}$ periods in the first reduction. The float conversion of that remainder has error below $10^{-15}$.
- Conversion and addition of the two limbs introduce error below $1.12\times10^{-12}$ in $\delta$. Multiplication by $\omega$ magnifies this by at most 728. Rounding $\omega$ contributes at most $3\times10^{-10}$ to $\omega\delta$. For the remaining operations, local phases stay below $3.65\times10^6$ and the rounded period differs from $2\pi$ by less than $6\times10^{-16}$. Fewer than $6\times10^5$ additional periods are subtracted. These estimates give phase error below $4\times10^{-9}$.
- The sine and cosine polynomials have degrees 25 and 24. On the checked local domain $|t|<3.142$, the alternating Taylor remainders are bounded by $3.142^{27}/27!$ and $3.142^{26}/26!$. The absolute sum of polynomial terms is less than 24 even for the sine polynomial after multiplying by $t$. Horner evaluation and rounded coefficients contribute less than $10^{-12}$. Including the phase error, both functions have absolute error below $5\times10^{-9}$.
- Write $v=g_0^2/(2\alpha)$ and $q=(\gamma^2-g_0^2)/(2\alpha)$. The scalar factor $e^{-v}$ uses the degree-17 Taylor polynomial at $-v/16$ followed by four squarings. Here $v\le 10.5$. The block correction $e^{-q}$ uses its cubic Taylor polynomial with $0\le q\le0.009$. Its omitted term has absolute size at most $0.009^4/24$. With coefficient and operation rounding, the combined absolute weight error is below $4\times10^{-10}$.
- For $c=1/\rho+1/(\omega\rho^2)$ with $\rho=1/2+i\gamma$, the real and imaginary coefficients are evaluated by rational expressions. Their rounding errors are bounded by $10^{-12}/\gamma$ from the displayed formulas in the implementation and the bounds $\gamma>14$ and $\omega>727$. Also $|\mathrm{Re}c|+|\mathrm{Im}c|\le\sqrt{2}(1+1/(\omega\gamma))/\gamma$. Including the bounded sine, cosine and weight errors, rounding the local sum with at most 50,000 terms, and the LMFDB ordinate radius of $2.5\times10^{-31}$, the error budget is $2\times10^{-8}/\gamma$ per conjugate pair.

Revers’s Lemma 4.8 bounds the sum of these allowances over the positive ordinates through the replay height. A further $10^{-14}$ covers the two levels of compensated summation. At the final cutoff the allowance is about $8.1\times10^{-7}$. The replay evaluates the margin from its own lower enclosure and uses Revers’s Gaussian tail bound at the actual cutoff. [Revers, Lemma 4.8 and Equation (5.5)](https://arxiv.org/html/2501.04488v2#S4)

## Validation

| Check | Observed result |
| --- | --- |
| Teal Sea control | $1.0028892744171265$ from 22 million zeros |
| Smith parameters through height $67{,}646{,}000$ | $1.0022437819281689$ from 163,559,197 zeros |
| Direct 60-digit check on the first archive block | Absolute difference below $1.65\times10^{-11}$ |
| Direct 60-digit checks on the first and last blocks of the final file | Absolute differences below $1.02\times10^{-23}$ |
| Final archive header | Height $30{,}610{,}046{,}000$, cumulative count $103{,}800{,}788{,}359$ |

The direct checks decoded the 13-byte offsets with Python integers and used the complex exponential formula in mpmath. The prefix through height $67{,}646{,}000$ took about 34 seconds to compute from 2.13 GB of cached data. The full run resumed from its SQLite checkpoint after a download failure. The completed checkpoint passed an integrity check, and its files, hashes, and heights matched the manifest sequence.

## Running the replay

Install NumPy and mpmath. Use the [public checksum manifest](https://beta.lmfdb.org/data/riemann-zeta-zeros/md5.txt) and a directory for the [binary files](https://www.lmfdb.org/zeros/zeta/). Start with the published 22-million-zero control from Teal Sea. The control uses the same decoder, arithmetic, and index checks at different Gaussian parameters. [Teal Sea, §7](https://zeta.teal-sea.com/read/hunts/skewes_enclosure/RESULTS.html)

```sh
python audit_sum.py data/platt data/platt-md5.txt data/audit-control.sqlite --control
python audit_sum.py data/platt data/platt-md5.txt data/audit-production.sqlite --through 67646000
python audit_sum.py data/platt data/platt-md5.txt data/audit-production.sqlite --download
```

Each source file is committed to SQLite only after its checksum, block structure, and contribution have been checked. The final command resumes from that checkpoint. It reads cached files if present and obtains missing files from LMFDB using curl. Each downloaded file uses a temporary path and is discarded after its contribution is committed. The complete archive transfers about 1.35 TB, while the temporary data occupy roughly one file at a time.
