"""Regression test: centroid_Calculation with a non-string cluster_col.

KNN2 always names its window columns as strings (it strips a
"{cluster_col}__" prefix off pandas get_dummies() column names). Before this
fix, centroid_Calculation indexed those string-named columns with the raw
`adata.obs[cluster_col].unique()` values, which raised a KeyError whenever
cluster_col held integers (a common encoding for cell-type ids).
"""

import anndata as ad
import numpy as np
import pandas as pd

from mingl.tl.centroids import centroid_Calculation


def _make_int_cluster_adata() -> ad.AnnData:
    rng = np.random.default_rng(0)
    n = 320  # KNN2's default ks include 300, so n must exceed that
    obs = pd.DataFrame(
        {
            "x": rng.uniform(0, 100, n),
            "y": rng.uniform(0, 100, n),
            "unique_region": ["R1"] * n,
            "cell_type": rng.integers(1, 4, n),  # integer labels: 1, 2, 3
            "neighborhood": rng.choice(["A", "B"], n),
        }
    )
    return ad.AnnData(X=np.zeros((n, 0), dtype=np.float32), obs=obs)


def test_centroid_calculation_accepts_integer_cluster_labels():
    adata = _make_int_cluster_adata()

    result = centroid_Calculation(
        adata, k=5, cluster_col="cell_type", neighborhood_col="neighborhood", region_col="unique_region"
    )

    assert set(result.obs_names) == {"A", "B"}
    # one mean/std column per distinct integer cell type
    expected_features = {f"{ct}_mean" for ct in (1, 2, 3)} | {f"{ct}_std" for ct in (1, 2, 3)}
    assert expected_features.issubset(set(result.var_names))


def test_centroid_calculation_string_and_int_labels_agree():
    adata = _make_int_cluster_adata()
    adata_str = adata.copy()
    adata_str.obs["cell_type"] = adata_str.obs["cell_type"].astype(str)

    result_int = centroid_Calculation(
        adata, k=5, cluster_col="cell_type", neighborhood_col="neighborhood", region_col="unique_region"
    )
    result_str = centroid_Calculation(
        adata_str, k=5, cluster_col="cell_type", neighborhood_col="neighborhood", region_col="unique_region"
    )

    np.testing.assert_allclose(
        result_int.to_df().reindex(index=result_str.obs_names, columns=result_str.var_names).values,
        result_str.to_df().values,
    )


def test_centroid_calculation_drops_cells_missing_neighborhood_label():
    """A NaN neighborhood_col value must never become a centroid name.

    Left in, it becomes a NaN-named row in the returned AnnData's obs_names,
    which downstream code that copies centroid names into new .obs columns
    (gpu_gmm_probability does this) cannot write to an h5ad file.
    """
    rng = np.random.default_rng(0)
    n = 320
    neighborhood = rng.choice(["A", "B"], n).astype(object)
    neighborhood[5] = np.nan  # one cell has no neighbourhood assignment
    obs = pd.DataFrame(
        {
            "x": rng.uniform(0, 100, n),
            "y": rng.uniform(0, 100, n),
            "unique_region": ["R1"] * n,
            "cell_type": rng.choice(["T", "B"], n),
            "neighborhood": neighborhood,
        }
    )
    adata = ad.AnnData(X=np.zeros((n, 0), dtype=np.float32), obs=obs)

    result = centroid_Calculation(
        adata, k=5, cluster_col="cell_type", neighborhood_col="neighborhood", region_col="unique_region"
    )

    assert set(result.obs_names) == {"A", "B"}
    assert not result.obs_names.isna().any()
