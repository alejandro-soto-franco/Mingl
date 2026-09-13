import anndata as ad
import numpy as np
import pandas as pd

from .knn2 import KNN2


def centroid_Calculation(
    adata: ad.AnnData,
    *,
    k: int = 10,
    cluster_col: str = "cell_type",
    neighborhood_col: str = "neighborhood",
    region_col: str = "unique_region",
    store_key: str | None = None,
) -> ad.AnnData:
    """
    Compute per-neighborhood mean and std of cell-type counts in k-NN windows.

    This is the same computation you had before, but:

    - keeps indices aligned (no reset_index)
    - allows choosing k
    - optionally stores the centroid AnnData in `adata.uns[store_key]`
      so you can reuse it for many plots.

    Parameters
    ----------
    adata
        AnnData with:
          - `.obs[cluster_col]` (cell type labels)
          - `.obs[neighborhood_col]` (neighborhood assignment)
          - spatial info needed for KNN (x/y & region columns).
    k
        Neighborhood size to use from the KNN windows (must be in the `ks`
        list passed to KNN; default 10).
    cluster_col
        Column in `adata.obs` with cell-type labels.
    neighborhood_col
        Column in `adata.obs` with neighborhood IDs.
    store_key
        If not None, store the centroid AnnData in `adata.uns[store_key]`.

    Returns
    -------
    AnnData
        AnnData with:
          - obs: neighborhoods
          - var: centroid features (means/stds per cell type)
          - X: numeric matrix (n_neighborhoods x n_features)

    Notes
    -----
    A cell missing `cluster_col` or `neighborhood_col` cannot contribute a
    labelled feature count or be assigned to a named centroid, so such rows
    are dropped before anything else runs (a real gap in practice: not every
    cell carries every level of a hierarchical annotation). Left in, a NaN
    `neighborhood_col` value becomes a NaN-named centroid, which downstream
    code that copies centroid names into new columns -- as
    :func:`mingl.tl.gmm_gpu.gpu_gmm_probability` does -- cannot write to an
    h5ad file.
    """
    n_before = adata.n_obs
    adata = adata[adata.obs[[cluster_col, neighborhood_col]].notna().all(axis=1)].copy()
    n_dropped = n_before - adata.n_obs
    if n_dropped:
        print(f"centroid_Calculation: dropped {n_dropped} cells missing {cluster_col!r}/{neighborhood_col!r}.")

    # get KNN windows
    windows = KNN2(adata, region_key=region_col, cluster_col=cluster_col)
    if k not in windows:
        raise ValueError(f"k={k} not in available ks from KNN: {list(windows.keys())}")

    windows_k = windows[k]

    # ensure we have the cluster column on windows_k
    windows_k[cluster_col] = adata.obs[cluster_col]

    # use obs directly; keep original indices (no reset_index)
    filtered_cells = adata.obs.copy()

    # cell types -> columns we created in KNN. KNN2's window columns are
    # always strings (it strips a "{cluster_col}__" prefix off pandas
    # get_dummies() column names), so cast here too: indexing windows_k with
    # the raw (e.g. integer) unique values would otherwise raise a KeyError
    # whenever cluster_col holds a non-string dtype.
    cell_type_columns = pd.Index(adata.obs[cluster_col].unique()).astype(str)
    windows_k[cell_type_columns] = windows_k[cell_type_columns].astype("float32")

    neighborhoods_to_loop = adata.obs[neighborhood_col].unique()
    all_results = []

    for neighborhood in neighborhoods_to_loop:
        # cells in this neighborhood (indices are original obs index)
        filtered_neighborhood_df = filtered_cells[filtered_cells[neighborhood_col] == neighborhood]
        cell_numbers_in_neighborhood = filtered_neighborhood_df.index.values

        # take matching rows from windows_k
        matching_cells_df = windows_k.loc[cell_numbers_in_neighborhood]

        mean_std_results = {neighborhood_col: neighborhood}

        for column in cell_type_columns:
            if column in matching_cells_df.columns:
                col_values = matching_cells_df[column]
                mean_std_results[f"{column}_mean"] = col_values.mean()
                mean_std_results[f"{column}_std"] = col_values.std()

        all_results.append(mean_std_results)

    # neighborhood × feature table
    results_df = pd.DataFrame(all_results).set_index(neighborhood_col)

    feature_cols = results_df.columns.tolist()
    X = results_df[feature_cols].to_numpy(dtype=np.float64)

    # obs: neighborhoods
    obs = pd.DataFrame(index=results_df.index)
    obs[neighborhood_col] = results_df.index

    # var: feature names
    var = pd.DataFrame(index=feature_cols)

    centroid_adata = ad.AnnData(X=X, obs=obs, var=var)

    if store_key is not None:
        adata.uns[store_key] = centroid_adata

    return centroid_adata
