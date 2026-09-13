"""CLI: reconstruct the notebooks' filtered intestine cell universe.

The deleted notebooks' own stored output reports 2,512,002 cells
(`fig2_intestine_tissue_unit.md`'s `n_obs`); the public Dryad file has
2,603,217. The notebooks read a further donor-metadata-merged, filtered
local copy whose exact derivation is not recoverable (see README's "Data
availability"). The closest reconstruction found: dropping the 91,032 cells
with no `Tissue Unit` label leaves 2,512,185 -- 183 cells (0.007%) away.

Applying that one filter once, here, and having every intestine rule read
its output instead of the raw file keeps fig2/fig3/fig4 on the same,
consistent, reconstructed cell universe, so their numbers are comparable to
each other and to the notebooks' own -- rather than each rule silently
using a different subset depending on which column it happens to need.
"""

from __future__ import annotations

import argparse


def main() -> None:
    """Parse CLI arguments and write the filtered intestine CSV."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-csv", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--filter-col", required=True, help="Column whose NaN rows are dropped.")
    args = parser.parse_args()

    import pandas as pd

    df = pd.read_csv(args.raw_csv)
    n_before = len(df)
    df = df.dropna(subset=[args.filter_col])
    n_dropped = n_before - len(df)
    print(f"Dropped {n_dropped} of {n_before} cells missing {args.filter_col!r}; {len(df)} remain.")
    df.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
