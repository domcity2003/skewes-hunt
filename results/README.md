# Numerical replay outputs

The [summary](zero-sum.json) records the run inputs, checksum manifest fingerprint, outward-rounded enclosure, and supported crossover interval. The [compressed CSV](file-sums.csv.gz) lists one row per source file. Each row records its archive name, published MD5 checksum, ending height, cumulative count, computed contribution, and byte count.

Run this from the repository root to reconstruct the total from the saved contributions. It needs only the Python standard library.

```python
import csv
import gzip
import math

with gzip.open("results/file-sums.csv.gz", "rt", encoding="ascii", newline="") as stream:
    rows = list(csv.DictReader(stream))
print(len(rows), rows[-1]["cumulative_zeros"])
print(math.fsum(float(row["score"]) for row in rows))
```

The expected lines are `14580 103800788359` and `1.0028636381592315`. This reconstructs the reported point estimate. The [full replay](../audit.md#running-the-replay) verifies each contribution from the original LMFDB files, and the [arithmetic notes](../audit.md#arithmetic) explain how to bound the accumulated numerical error.
