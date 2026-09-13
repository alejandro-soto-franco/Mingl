"""CLI: compute and cache the k=10 cell-type composition window once.

fig6's n-cluster sweep (`run_mingl_over_n_clusters`) re-clusters the SAME
per-cell composition features at every candidate cluster count; it does not
recompute a KNN2 window per n. Computing that window once here, and having
every per-n job (`run_n_neighborhood_cluster.py`) read this cached output
instead of recomputing it, is what makes staging the sweep over 30-minute,
per-computation jobs meaningful: the expensive, real-data-scale step (a
k=10 KNN2 pass over every surviving cell) happens exactly once.
"""

from __future__ import annotations

import argparse

import anndata as ad
import numpy as np
import pandas as pd

from mingl.tl.knn2 import KNN2


def main() -> None:
    """Parse CLI arguments and cache the k-window composition matrix."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-csv", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--cluster-col", required=True)
    parser.add_argument("--region-key", required=True)
    parser.add_argument("--x-key", default="x")
    parser.add_argument("--y-key", default="y")
    parser.add_argument("--k", type=int, required=True)
    args = parser.parse_args()

    usecols = [args.x_key, args.y_key, args.region_key, args.cluster_col]
    dtype = {args.region_key: "category", args.cluster_col: "category"}
    df = pd.read_csv(args.raw_csv, usecols=usecols, dtype=dtype)
    df[args.x_key] = df[args.x_key].astype(np.float32)
    df[args.y_key] = df[args.y_key].astype(np.float32)
    cells = ad.AnnData(X=np.zeros((len(df), 0), dtype=np.float32), obs=df)

    windows = KNN2(
        cells,
        x_key=args.x_key,
        y_key=args.y_key,
        region_key=args.region_key,
        cluster_col=args.cluster_col,
        ks=(args.k,),
        keep_obs_cols=[args.x_key, args.y_key, args.region_key],
    )
    window_df = windows[args.k]
    # KNN2 always prepends x_key/y_key/region_key/cluster_col ahead of the
    # dummy-encoded feature sums (cluster_col along with them, regardless of
    # keep_obs_cols): exclude all four, not just the three passed above.
    non_feature_cols = {args.x_key, args.y_key, args.region_key, args.cluster_col}
    feature_cols = [c for c in window_df.columns if c not in non_feature_cols]

    out = ad.AnnData(
        X=np.zeros((len(window_df), 0), dtype=np.float32),
        obs=window_df[[args.x_key, args.y_key, args.region_key]].copy(),
    )
    out.obsm[f"knn_windows_k{args.k}"] = window_df[feature_cols].astype("float32").to_numpy()
    out.uns["knn_windows"] = {"k": args.k, "cols": list(feature_cols), "cluster_col": args.cluster_col}
    out.write_h5ad(args.output)
    print(f"Cached k={args.k} windows for {out.n_obs} cells, {len(feature_cols)} feature columns.")


if __name__ == "__main__":
    main()
