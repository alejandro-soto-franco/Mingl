import time
from collections.abc import Sequence

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors


def _sum_windows(window: np.ndarray, dists_k: np.ndarray, max_distance: float | None) -> np.ndarray:
    """Sum a ``(n_cells, k, n_features)`` window array over the neighbour axis.

    This is the one place both :func:`KNN2` and
    :meth:`mingl.tl.grad.Neighborhoods.k_windows` reduce a k-NN window to a
    per-cell feature sum, so a distance threshold behaves identically for
    both. When ``max_distance`` is ``None`` every neighbour counts (the
    original, unthresholded behaviour). Otherwise, neighbours farther than
    ``max_distance`` (per ``dists_k``, matched neighbour-for-neighbour with
    ``window``'s second axis) are excluded from the sum for that cell.
    """
    if max_distance is None:
        return window.sum(axis=1)
    mask = dists_k > max_distance
    mask3 = np.repeat(mask[:, :, np.newaxis], window.shape[2], axis=2)
    return np.ma.array(window, mask=mask3).sum(axis=1).data


def KNN2(
    adata: ad.AnnData,
    *,
    x_key: str = "x",
    y_key: str = "y",
    region_key: str = "unique_region",
    cluster_col: str = "cell_type",
    ks: Sequence[int] = (5, 10, 20, 100, 300),
    keep_obs_cols: Sequence[str] | None = None,  # <-- NEW
    max_distance: float | None = None,
) -> dict[int, pd.DataFrame]:
    """
    Compute cell-type neighborhood windows using a k-NN in (x,y) per region.

    NEW:
    ----
    keep_obs_cols:
        Optional list/tuple of columns from adata.obs to carry through into each
        output window dataframe (in addition to region_key and cluster_col).
    max_distance:
        Optional maximum neighbour distance. When given, neighbours farther
        than this are excluded from a cell's window sum instead of counted
        (the same thresholding as ``Neighborhoods.k_windows`` in
        :mod:`mingl.tl.grad`, implemented via the shared ``_sum_windows``
        helper). ``None`` (the default) reproduces the original, unthresholded
        behaviour exactly.

    Returns
    -------
    dict[int, pd.DataFrame]
        Mapping k -> DataFrame of neighborhood summaries.
    """
    ks = list(ks)
    n_neighbors = max(ks)

    adata = adata.copy()
    if adata.obs.columns.has_duplicates:
        counts = {}
        new_cols = []
        for c in adata.obs.columns:
            if c not in counts:
                counts[c] = 0
                new_cols.append(c)
            else:
                counts[c] += 1
                new_cols.append(f"{c}__dup{counts[c]}")
        adata.obs.columns = pd.Index(new_cols)

    # Default keep set: keep x/y if present (without changing your logic elsewhere)
    if keep_obs_cols is None:
        keep_obs_cols = []

    # Build the final "kept" columns list (preserve user order, de-dup)
    keep_cols: list[str] = []
    for c in list(keep_obs_cols) + [x_key, y_key, region_key, cluster_col]:
        if c not in keep_cols:
            keep_cols.append(c)

    # Filter to only those that exist in obs (avoid KeyError if user passes extras)
    keep_cols_existing = [c for c in keep_cols if c in adata.obs.columns]

    # Dummy encode cell types (same logic)
    # Dummy encode cell types
    dummies = pd.get_dummies(adata.obs[cluster_col]).add_prefix(f"{cluster_col}__")
    adata.obs = pd.concat([adata.obs, dummies], axis=1)
    print("Example dummy cols:", list(dummies.columns[:5]))

    # Use the dummy columns as sum_cols (these ARE real columns in obs)
    sum_cols = dummies.columns.to_numpy()
    values = adata.obs[sum_cols].to_numpy()

    tissue_group = adata.obs.groupby(region_key)
    exps = list(adata.obs[region_key].unique())

    def get_windows(job, n_neighbors: int) -> tuple[np.ndarray, np.ndarray]:
        _start_time, _idx, tissue_name, indices = job

        region = adata.obs[region_key]
        # If region_key hits duplicate columns, pandas returns a DataFrame — take the first.
        if isinstance(region, pd.DataFrame):
            region = region.iloc[:, 0]

        mask = region.to_numpy() == tissue_name
        tissue = adata[mask].copy()  # avoid AnnData view machinery

        coords = tissue.obs[[x_key, y_key]].values
        fit = NearestNeighbors(n_neighbors=n_neighbors).fit(coords)
        distances, neighbor_idx = fit.kneighbors(coords)

        args = distances.argsort(axis=1)
        add = np.arange(neighbor_idx.shape[0]) * neighbor_idx.shape[1]
        sorted_indices = neighbor_idx.flatten()[args + add[:, None]]
        sorted_dists = distances.flatten()[args + add[:, None]]

        obs_index = np.asarray(tissue.obs.index.to_numpy())
        neighbors = obs_index[sorted_indices]
        return neighbors.astype(np.int32), sorted_dists

    tissue_chunks = [
        (time.time(), exps.index(tissue_name), tissue_name, indices)
        for tissue_name, indices in tissue_group.groups.items()
        for indices in np.array_split(indices, 1)
    ]

    tissues = [get_windows(job, n_neighbors) for job in tissue_chunks]

    out_dict = {}
    for k in ks:
        for (neighbors, dists), job in zip(tissues, tissue_chunks):
            tissue_name = job[2]
            indices = job[3]
            chunk = np.arange(len(neighbors))

            window = values[neighbors[chunk, :k].flatten()]
            window = window.reshape(len(chunk), k, len(sum_cols))
            summed = _sum_windows(window, dists[:, :k], max_distance)
            out_dict[(tissue_name, k)] = (summed.astype(np.float16), indices)

    windows: dict[int, pd.DataFrame] = {}
    for k in ks:
        dfs = []
        for exp in exps:
            data, idx = out_dict[(exp, k)]
            # Keep your original behavior of indexing by idx.astype(int)
            pretty_cols = [c.removeprefix(f"{cluster_col}__") for c in sum_cols]
            df_k = pd.DataFrame(data, index=idx.astype(int), columns=pretty_cols)

            dfs.append(df_k)

        window_df = pd.concat(dfs, axis=0)

        # Ensure the index types match (your existing guard)
        if window_df.index.dtype != adata.obs.index.dtype:
            window_df.index = window_df.index.astype(str)

        window_df = window_df.loc[adata.obs.index]

        # NEW: prepend kept columns (including x/y if requested / available)
        window_df = pd.concat(
            [adata.obs[keep_cols_existing], window_df],
            axis=1,
        )

        windows[k] = window_df

    return windows
