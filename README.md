# Skewes crossover search

Started November 20, 2026 (US Central) with Codex / GPT-6. The target is a proof that the actual prime count exceeds the logarithmic integral somewhere below Smith's thesis bound, `log(x) = 727.9513354023`. The logarithm is natural except in the Table 5 scan.

The [finite-height calculation](unconditional-bound.md) places a crossover at or below `log(x) = 727.95133540685` using Smith and Platt's published zero sum. The arithmetic can be replayed locally. The source zero sum and literature priority still need independent review.

## Next experiment

Smith's Section 10 reports a signal near `727.951332982` using about 100 billion zeros, with an error allowance larger than the signal. The thesis does not report the parameters or raw margin for that run. A search of the public thesis deposit and indexed code did not turn up a run log. Those values would tell us how much of the error we need to remove and which term to target. [Smith, Section 10](https://etheses.whiterose.ac.uk/16409/1/The%20hunt%20for%20Skewes%27%20number.pdf)

The [resumable sweep](stream_sweep.py) reads consecutive LMFDB files, verifies their checksums and block continuity, and saves the accumulated sums after each file. It checks the positive control at exactly 22 million zeros, inside a file. At each cutoff it records the weighted zero sum and each of Revers's error terms, including the omitted-zero allowance `R4`. This uses Revers's prime-side bound at every center. A candidate with a positive exploratory margin would need interval arithmetic on the original Platt ordinates. [LMFDB data](https://www.lmfdb.org/zeros/zeta/) · [Revers's explicit bounds](https://arxiv.org/html/2501.04488v2#S3)

The file names and MD5 digests used in the experiment are recorded in [stream-sweep.json](stream-sweep.json). Download the files from the LMFDB bulk archive into `data/platt` and save the [manifest](https://beta.lmfdb.org/data/riemann-zeta-zeros/md5.txt) as `data/platt-md5.txt`. Use the ending height of the last downloaded file for `--through`. The checkpoint stays in `data/`. Rerun the command with a higher ending height after adding the next consecutive files. One OpenBLAS thread ran faster on the Windows machine used here.

```powershell
$env:OPENBLAS_NUM_THREADS = '1'
python stream_sweep.py data/platt data/platt-md5.txt data/stream-checkpoint.json stream-sweep.json --through 65546000
```

The November 21, 2026 run (Codex / GPT-6) reached 158,152,688 zeros in 35 files, testing 25 centers and 50 Gaussian widths. The 22-million-zero control returned `1.002889274491` with a positive full-RH margin. The highest raw sum across all early-target checkpoints was `1.002313399511`, below the prime-side threshold of about `1.002762685`. At the final cutoff, the best finite-height margin was `-0.0008806563` at the named target. The omitted-zero allowance `R4` in that row was `0.0000848369`, smaller than the raw deficit of about `0.000794271`. The scan used float64 arithmetic and leaves the earlier crossover unresolved.

An independent evaluation using the explicit sine and cosine formula agreed at the 99,480,453-zero checkpoint within `1.4e-10`. The decoder matched the previous 23,293,070 ordinates bit for bit, and replaying the completed checkpoint reproduced the report byte for byte.

At the early target, reusing the smoothing `alpha = 4.5e19` and cutoff `T = 30610046000` from Smith's confirmed interval gives a screening threshold of about `1.002772957` for the weighted sum. This uses the largest trial window allowed by the endpoint target and Revers's finite-height bound with `A = 2e11`. The earlier run's smoothing, cutoff, and sum remain unknown. [Smith, Section 9](https://etheses.whiterose.ac.uk/id/eprint/16409/) · [Revers, Theorem 3.1](https://arxiv.org/html/2501.04488v2#S3)

## Reproduce the first pass

Use Python 3.12, NumPy, and mpmath. `data/` is a local cache. The setup below fetches the same source files used for the run.

```powershell
python -m pip install numpy mpmath
New-Item -ItemType Directory -Path data -Force
curl.exe -fL https://raw.githubusercontent.com/mesito/Speiser_Lehmer_identity/main/data/zeros6.gz -o data/zeros6.gz
git clone --filter=blob:none --sparse https://github.com/teal-sea/zeta-lab.git data/zeta-lab
git -C data/zeta-lab checkout a401821c9608526425f3d4cf7ff8e59bc37e1c79
git -C data/zeta-lab sparse-checkout set hunts/skewes_enclosure
python search.py data/zeros6.gz --early --output pilot.json
python replay_grid.py data/zeta-lab/hunts/skewes_enclosure grid-replay.json
python replay_smith.py
```

The compressed Odlyzko input has SHA-256 `6acacb4707c429bf368f32823a8eef0e0737b6318d3e144c78975a4986d155da`. The benchmark and Table 5 search use exploratory float64 arithmetic. The published cutoff separates the first two million zeros from the rest of the table. The benchmark matches the reference sum within `4.64e-8` and reproduces the switch in the sign of the estimated margin as the interval tightens.

## Bounding the spaces between grid points

`replay_grid.py` adds Teal Sea's recorded Arb enclosures using exact rationals. It checks shard indices, inter-shard ordinate order, cutoff placement, and agreement with published totals. It relies on the upstream isolation and Turing count. The column of estimated margins uses floating point. The exclusion uses exact rational inequalities and applies specifically to Revers's stated lower bound with `alpha = 6e12` and `T = 10379600`.

The margin is at most `-0.00108134236` throughout the center range from `727.951333792` to `727.9513354023`, for every admissible positive `eta` at those fixed values of `alpha` and `T`. A finer grid or a different `eta` in this range cannot yield the desired positive lower bound. The next search should vary the Gaussian width and the zero cutoff together. The [Saouter–Trudgian–Demichel weighted kernel](https://arxiv.org/html/2501.04488v2#S7) is another way to change that tradeoff.

For positive `eta`, Revers's first error term alone is larger than `2/omega`. If the finite zero sum `S(omega)` is too high for `-1 - S(omega) - 2/omega` to be positive, changing `eta` at that center cannot give a positive lower bound by this formula. This criterion also applies on the full-RH branch where `R6` is omitted.

The grid spacing is `5e-8`. For each segment, four adjacent grid values determine an interpolating cubic. Its lowest Bernstein coefficient bounds it from below on that segment. The absolute interpolation remainder is bounded by `M * step^4 / 24`, since the product of distances to four equally spaced nodes is at most `step^4` between the outer nodes. In units of the step, this follows from `(z - 9/4)(z - 1/4)` with `z` in `[0, 9/4]`.

For `omega > 700`, a bound `M = 9e26` suffices for the fourth derivative of the finite sum. The fourth derivative of the `1/rho` part has magnitude at most `2N * max(g^3 exp(-g^2/(2*alpha)))`, where `N = 22,000,000`. The maximum occurs at `g = sqrt(3*alpha) < 4.25e6` and is less than `(4.25e6)^3/4` since `exp(3/2) > 4`. For the `1/(omega*rho^2)` part, the product rule adds at most this same bound times `1/700 + 4/700^2 + 12/700^3 + 24/700^4 + 24/700^5 < 0.002`. This uses `g > 1` and `abs(rho) >= g`. Together the bounds are less than `M`.

## Candidate improvement from Smith's tabulated sum

`replay_smith.py` holds Smith's center `727.951335401`, `alpha = 4.5e19`, and `M = T = 30610046000` fixed. It uses Smith's Theorem 4.1, which assumes RH. The two sums in equation (26) are independent of `eta`. Appendix C, Table 4 reports their combined contribution as `1.00288049783641096204`. The replay uses the conservative lower value `1.0028804978364`. Section 9 gives the original value of `C` at `eta = 1.3e-9` as `-1.002762659…`. The integrand in (27) is positive over these intervals before the leading minus sign, so reducing `eta` makes `C` less negative and the lower bound `C > -1.002763` still applies.

At `eta = 7.5e-10`, interval arithmetic gives total error below `0.000080` and a lower bound on `I` above `0.000037`. Thus the theorem and reported inputs yield a candidate interval `[727.95133540025, 727.95133540175]` in `log(x)` coordinates. Its upper endpoint is `5.5e-10` earlier in `log(x)` than Theorem 9.1's endpoint. This calculation re-evaluates the error terms. It uses Smith and Platt's reports as inputs for the zero sum and original integral value. A separate audit of those source computations and priority remains open. Smith also reports an earlier signal near `727.951332982` that was smaller than his error allowance. That is a specific target for a larger gain.

## Follow-up using the original Platt files

The first nine Platt files, ending with `zeros_8846000.dat`, contain 23,293,070 consecutive ordinates. The original files total 303 MB. `platt.py` checks the official MD5 manifest, continuity of the block indices and heights, and the order and range of the decoded ordinates. Its NumPy output rounds ordinates to float64 for scouting. The Gaussian sweep independently reproduces Teal Sea's 22-million-zero control and positive benchmark margin. At the earlier candidate, its best estimated margin is `-0.0020551`. Even the largest sampled zero sum, `1.0015162`, is below the prime-side allowance near `1.002763` before the analytic errors are deducted.

`compact.py` tries the normalized nonnegative weight `(1-(y/eta)^2)^p` supported on `[-eta,eta]`. The transform is `(2p+1)!! j_p(z)/z^p` where `z = eta*gamma`. It uses Revers (4.10)–(4.12) for the prime-side and expansion allowances and bounds the omitted zero sum using the finite Hankel expansion and Saouter–Trudgian–Demichel Lemma 2.7. It screened eight orders at 13 centers. Its largest estimated margin was `-0.002224`. These sweeps leave the earlier crossover unresolved. [Bessel integral](https://dlmf.nist.gov/10.54.E1) · [Hankel expansion](https://dlmf.nist.gov/10.49) · [Revers](https://arxiv.org/html/2501.04488v2) · [Saouter–Trudgian–Demichel](https://openresearch-repository.anu.edu.au/server/api/core/bitstreams/3bbcd8d2-860d-434a-8708-eb7c12e33986/content)

For locally downloaded files in `data/platt` and the official `md5.txt` at `data/platt-md5.txt`, run the following commands. The compact-weight sweep also needs SciPy.

```powershell
python platt.py data/platt data/platt-md5.txt data/platt.npy
python earlier.py data/platt.npy earlier-sweep.json
python compact.py data/platt.npy compact-sweep.json
```

Separately, [the finite-height bound](unconditional-bound.md) retains all six errors of Revers's theorem and uses a published verification only up to a finite height. Its upper endpoint is `727.95133540685`. Run `python unconditional.py` to reproduce both the interval calculation and the exact rational cross-check.

## Sources

- [Revers, Theorem 3.1 and Table 5](https://arxiv.org/html/2501.04488v2)
- [Smith, The hunt for Skewes’ number, Theorem 9.1](https://etheses.whiterose.ac.uk/id/eprint/16409/)
- [Teal Sea’s source code](https://github.com/teal-sea/zeta-lab/tree/main/hunts/skewes_enclosure) and [computation report](https://zeta.teal-sea.com/read/hunts/skewes_enclosure/RESULTS.html)
- [Odlyzko’s original zero tables](https://www-users.cse.umn.edu/~odlyzko/zeta_tables/)
- [LMFDB’s certified Platt dataset](https://www.lmfdb.org/zeros/zeta/) and [binary reader](https://github.com/LMFDB/lmfdb/blob/main/lmfdb/zeros/zeta/platt_zeros.py)
