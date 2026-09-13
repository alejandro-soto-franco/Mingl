"""Regression tests for KNN2's optional distance threshold.

These cover the consolidation of KNN2 (mingl.tl.knn2) and
Neighborhoods.k_windows (mingl.tl.grad) onto one shared summation helper,
``mingl.tl.knn2._sum_windows``: both now thin-wrap it, and a distance
threshold behaves identically through either entry point.
"""

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from mingl.tl.knn2 import KNN2, _sum_windows


def _make_line_adata() -> ad.AnnData:
    """Five cells on a line, one per unit distance, one region, two types."""
    obs = pd.DataFrame(
        {
            "x": [0.0, 1.0, 2.0, 3.0, 4.0],
            "y": [0.0, 0.0, 0.0, 0.0, 0.0],
            "unique_region": ["R1"] * 5,
            "cell_type": ["A", "B", "A", "B", "A"],
        },
        index=[str(i) for i in range(5)],
    )
    return ad.AnnData(X=np.zeros((5, 0), dtype=np.float32), obs=obs)


def test_sum_windows_none_matches_plain_sum():
    window = np.arange(24, dtype=np.float64).reshape(2, 3, 4)
    dists = np.array([[0.1, 0.2, 0.3], [0.1, 0.2, 0.3]])
    result = _sum_windows(window, dists, max_distance=None)
    np.testing.assert_allclose(result, window.sum(axis=1))


def test_sum_windows_excludes_far_neighbours():
    window = np.ones((1, 3, 2), dtype=np.float64)
    dists = np.array([[0.5, 1.5, 2.5]])
    result = _sum_windows(window, dists, max_distance=1.0)
    # only the first neighbour (distance 0.5) survives the threshold
    np.testing.assert_allclose(result, np.array([[1.0, 1.0]]))


def test_knn2_default_max_distance_is_unthresholded():
    adata = _make_line_adata()
    windows = KNN2(adata, ks=(2,))
    window_df = windows[2]
    # cell "0" (x=0): 2 nearest neighbours are cells "1" (B) and "2" (A)
    assert window_df.loc["0", "A"] == 1
    assert window_df.loc["0", "B"] == 1


def test_knn2_max_distance_excludes_far_neighbours():
    adata = _make_line_adata()
    unthresholded = KNN2(adata, ks=(2,), max_distance=None)[2]
    thresholded = KNN2(adata, ks=(2,), max_distance=0.5)[2]

    # KNN2's k-window for a cell includes the cell itself (distance 0);
    # cell "0"'s other k=2 neighbour ("1", type B) is at distance 1.0, so a
    # 0.5 cutoff excludes it and only the self-count (type A) remains.
    assert unthresholded.loc["0", ["A", "B"]].sum() == 2
    assert thresholded.loc["0", ["A", "B"]].sum() == 1
    assert thresholded.loc["0", "A"] == 1
    assert thresholded.loc["0", "B"] == 0


def test_knn2_raises_on_missing_columns():
    obs = pd.DataFrame({"x": [0.0], "y": [0.0], "unique_region": ["R1"]}, index=["0"])
    adata = ad.AnnData(X=np.zeros((1, 0), dtype=np.float32), obs=obs)
    with pytest.raises(KeyError):
        KNN2(adata, cluster_col="cell_type")
