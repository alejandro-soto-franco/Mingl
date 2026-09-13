"""CLI: generic k-window neighbourhood/tissue-unit/community pass.

The manuscript's intestine (fig2, fig4) and melanoma (fig2) figures share one
computation once their raw table is loaded: KNN2 windows -> per-neighbourhood
centroids -> Gaussian-mixture membership probabilities -> threshold count.
`centroid_Calculation` already calls KNN2 internally, so this script only
calls it and `cpu_gmm_probability`, parameterized entirely from config.yaml
(`datasets.<name>.neighborhood` / `.tissue_unit` / `.community`) rather than
the manuscript-specific defaults the library used to hardcode.

A plain CLI (not a Snakemake `script:`) so the calling rule can wrap it in
`flock` around the shared raw-data file: `--raw-csv` is read with `usecols`
restricted to only the columns this run needs, and cluster/neighbourhood/
region columns loaded as `category` dtype, to keep peak RSS well under the
~12 GB this machine's shared 30 GB is budgeted at.
"""

from __future__ import annotations

import argparse

import anndata as ad
import numpy as np
import pandas as pd

from mingl.tl.centroids import centroid_Calculation
from mingl.tl.edges import findPositives
from mingl.tl.gmm import cpu_gmm_probability


def _load_cells(
    raw_csv: str, *, x_key: str, y_key: str, region_key: str, cluster_col: str, neighborhood_col: str
) -> ad.AnnData:
    usecols = [x_key, y_key, region_key, cluster_col, neighborhood_col]
    dtype = {region_key: "category", cluster_col: "category", neighborhood_col: "category"}
    df = pd.read_csv(raw_csv, usecols=usecols, dtype=dtype)

    # A cell missing its cluster or neighbourhood label cannot contribute a
    # centroid or be scored against one; left in, a NaN neighbourhood name
    # becomes a NaN column name when probabilities are copied into .obs
    # below, which anndata's h5ad writer rejects outright. Drop them instead
    # of silently coercing NaN to a string category.
    n_before = len(df)
    df = df.dropna(subset=[cluster_col, neighborhood_col, region_key])
    n_dropped = n_before - len(df)
    if n_dropped:
        print(f"Dropped {n_dropped} of {n_before} cells missing {cluster_col!r}/{neighborhood_col!r}/{region_key!r}.")
        # KNN2 uses obs.index labels as positions into a plain numpy array
        # built from the same obs (values = adata.obs[sum_cols].to_numpy()),
        # so it silently assumes a contiguous 0..n-1 index. dropna() keeps
        # the original (now non-contiguous, and no longer 0..n-1) labels,
        # which raised "index ... out of bounds" once a surviving row's
        # original position exceeded the post-drop row count.
        df = df.reset_index(drop=True)

    # KNN2 groups by region_key with pandas' default observed=False, which
    # yields an (empty) group for every category the dtype still lists even
    # after the dropna above -- but adata.obs[region_key].unique() only
    # returns categories with actual rows. That mismatch made KNN2 raise
    # "<region> is not in list" for any region left with zero cells. Drop the
    # now-unused categories so both agree on which regions exist.
    for col in (region_key, cluster_col, neighborhood_col):
        df[col] = df[col].cat.remove_unused_categories()

    df[x_key] = df[x_key].astype(np.float32)
    df[y_key] = df[y_key].astype(np.float32)
    empty_x = np.zeros((df.shape[0], 0), dtype=np.float32)
    return ad.AnnData(X=empty_x, obs=df)


def main() -> None:
    """Parse CLI arguments and run the k-window neighbourhood-probability pipeline."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-csv", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--cluster-col", required=True)
    parser.add_argument("--neighborhood-col", required=True)
    parser.add_argument("--region-key", required=True)
    parser.add_argument("--x-key", default="x")
    parser.add_argument("--y-key", default="y")
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--ks", type=int, nargs="+", required=True)
    parser.add_argument("--threshold", type=float, default=0.25)
    args = parser.parse_args()

    cells = _load_cells(
        args.raw_csv,
        x_key=args.x_key,
        y_key=args.y_key,
        region_key=args.region_key,
        cluster_col=args.cluster_col,
        neighborhood_col=args.neighborhood_col,
    )

    centroids = centroid_Calculation(
        cells.copy(),
        k=args.k,
        cluster_col=args.cluster_col,
        neighborhood_col=args.neighborhood_col,
        region_col=args.region_key,
    )
    result = cpu_gmm_probability(
        cells,
        centroids,
        cluster_col=args.cluster_col,
        neighborhood_col=args.neighborhood_col,
        region_key=args.region_key,
        ks=tuple(args.ks),
        k=args.k,
    )
    findPositives(result, prob_key="neighborhood_probabilities", threshold=args.threshold)

    # Also copy the probability matrix into named .obs columns (matching
    # gpu_gmm_probability's own convention), since the network-graph rules
    # read per-neighbourhood probability columns from .obs, not .obsm.
    names = list(result.uns["neighborhood_probability_neighborhoods"])
    result.obs[names] = result.obsm["neighborhood_probabilities"]

    result.write_h5ad(args.output)


if __name__ == "__main__":
    main()
