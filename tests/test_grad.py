"""Regression test: Neighborhoods.k_windows with a non-contiguous index.

Same root cause as KNN2's fix (see tests/test_knn2.py): make_windows mapped
local, positional k-NN results to `tissue.index` labels, then used those
labels to index `values` -- a plain, positional numpy array. That only
gives the right answer when the DataFrame's index happens to already equal
0..n-1 in row order.
"""

import numpy as np
import pandas as pd

from mingl.tl.grad import Neighborhoods


def _make_frame(n: int) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "X:X": np.arange(n, dtype=np.float64),
            "Y:Y": np.zeros(n),
            "Exp": ["R1"] * n,
            "cell_type": ["A" if i % 2 == 0 else "B" for i in range(n)],
        }
    )


def _run(cells: pd.DataFrame, k: int = 3) -> pd.DataFrame:
    neigh = Neighborhoods(
        cells,
        ks=[k],
        cluster_col="cell_type",
        sum_cols=["A", "B"],  # add_dummies() one-hot-encodes cell_type into columns named by its own values
        keep_cols=["X:X", "Y:Y", "Exp"],
        neigh="Exp",
    )
    return neigh.k_windows()[k]


def test_neighborhoods_non_contiguous_index_matches_reset_index():
    full = _make_frame(12)
    filtered = full.drop(index=[2, 5, 8])  # index is now non-contiguous: 0,1,3,4,6,7,9,10,11
    reindexed = filtered.reset_index(drop=True)

    windows_non_contiguous = _run(filtered.copy())
    windows_reset = _run(reindexed.copy())

    dummy_cols = [c for c in windows_non_contiguous.columns if c not in ("X:X", "Y:Y", "Exp")]
    np.testing.assert_array_equal(
        windows_non_contiguous[dummy_cols].to_numpy(),
        windows_reset[dummy_cols].to_numpy(),
    )
