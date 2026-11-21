"""Decode consecutive Platt/LMFDB files for floating-point searches."""

import argparse
import hashlib
from pathlib import Path
import struct

import numpy as np


# Layout documented in Jonathan Bober's LMFDB platt_zeros.py.
HEADER = struct.Struct("<ddQQ")
GAP = np.dtype([("low", "<u8"), ("mid", "<u4"), ("high", "u1")])


def read_exact(stream, size):
    data = stream.read(size)
    if len(data) != size:
        raise ValueError(f"Truncated file: {stream.name}")
    return data


def ordinate_blocks(directory, manifest):
    """Yield contiguous blocks whose source files match the official manifest."""
    hashes = dict((name.lstrip("*"), digest) for digest, name in
                  (line.split() for line in manifest.read_text(encoding="ascii").splitlines() if line.strip()))
    files = sorted(directory.glob("zeros_*.dat"), key=lambda path: int(path.stem[6:]))
    if not files:
        raise ValueError("No ordinate files found")
    height, index = 14, 0
    for path in files:
        with path.open("rb") as stream:
            md5 = hashlib.file_digest(stream, lambda: hashlib.md5(usedforsecurity=False)).hexdigest()
            if md5 != hashes.get(path.name):
                raise ValueError(f"Checksum mismatch: {path.name}")
            stream.seek(0)
            size, = struct.unpack("<Q", read_exact(stream, 8))
            for _ in range(size):
                start, end, first, last = HEADER.unpack(read_exact(stream, HEADER.size))
                if (start, first) != (height, index) or not start < end or not first < last:
                    raise ValueError(f"Discontinuous or invalid block in {path.name}")
                gaps = np.frombuffer(read_exact(stream, GAP.itemsize * (last - first)), dtype=GAP)
                # Add the 104-bit gaps using two limbs. Convert only the exact
                # cumulative offsets to float64, after accounting for carries.
                low = np.cumsum(gaps["low"], dtype=np.uint64)
                upper = gaps["mid"].astype(np.uint64) + (gaps["high"].astype(np.uint64) << 32)
                upper[1:] += low[1:] < low[:-1]
                high = np.cumsum(upper, dtype=np.uint64)
                zeros = start + np.ldexp(high.astype(np.float64), -37) + np.ldexp(low.astype(np.float64), -101)
                if not start < zeros[0] or not zeros[-1] < end or not np.all(np.diff(zeros) > 0):
                    raise ValueError(f"Ordinates are outside the ordered block in {path.name}")
                yield zeros
                height, index = end, last


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    zeros = np.concatenate(list(ordinate_blocks(args.directory, args.manifest)))
    np.save(args.output, zeros, allow_pickle=False)
    print(f"Saved {len(zeros):,} consecutive ordinates, ending at height {zeros[-1]:.8f}")


if __name__ == "__main__":
    main()
