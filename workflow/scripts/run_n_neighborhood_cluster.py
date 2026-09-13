"""CLI: run fig6's n-cluster sweep for one n, against the cached windows.

`run_mingl_over_n_clusters` accepts any iterable for `n_range`; passing a
single-element list runs exactly one cluster count's KMeans fit plus
batched Gaussian-likelihood evaluation, identical to that same n inside the
full 1..50 loop (MiniBatchKMeans is seeded with random_state=0, so each n is
independent of loop order). One Snakemake job per n, each within the
30-minute-per-computation limit, is how the full sweep is staged; a
resumed `pixi run snakemake` skips whichever n{...}.json files already
exist.
"""

from __future__ import annotations

import argparse
import json

import anndata as ad

from mingl.tl.n_neighbors import run_mingl_over_n_clusters


def main() -> None:
    """Parse CLI arguments and run one n's clustering pass."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows", required=True, help="Path to compute_n_neighborhood_windows.py's output.")
    parser.add_argument("--output", required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--x-key", default="x")
    parser.add_argument("--y-key", default="y")
    parser.add_argument("--region-key", default="unique_region")
    args = parser.parse_args()

    adata = ad.read_h5ad(args.windows)
    meta = adata.uns["knn_windows"]
    feature_cols = list(meta["cols"])
    k = int(meta["k"])
    adata.obs[feature_cols] = adata.obsm[f"knn_windows_k{k}"]

    cluster_col = f"Neighborhood_{args.n}"
    summary_df = run_mingl_over_n_clusters(
        adata,
        knn_feature_cols=feature_cols,
        n_range=[args.n],
        output_col_template="Neighborhood_{}",
        return_per_cell=False,
        plot_summary=False,
        x_key=args.x_key,
        y_key=args.y_key,
        region_key=args.region_key,
    )

    row = summary_df.iloc[0]
    cluster_sizes = adata.obs[cluster_col].value_counts().sort_index()
    result = {
        "n": args.n,
        "avg_log_likelihood": float(row["avg_log_likelihood"]),
        "avg_assigned_probability": float(row["avg_assigned_probability"]),
        "n_cells": int(adata.n_obs),
        "cluster_sizes": {str(label): int(count) for label, count in cluster_sizes.items()},
    }
    with open(args.output, "w") as fh:
        json.dump(result, fh, indent=2)
    print(
        f"n={args.n}: avg_log_likelihood={row['avg_log_likelihood']:.4f}, "
        f"avg_assigned_probability={row['avg_assigned_probability']:.4f}"
    )


if __name__ == "__main__":
    main()
