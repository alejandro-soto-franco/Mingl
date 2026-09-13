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
    """Read only the columns this run needs, as category dtype where useful.

    Missing-label handling (a cell with no cluster/neighbourhood/region
    value) and non-contiguous-index safety are `centroid_Calculation`'s and
    `KNN2`'s own concerns now (see src/mingl/tl/centroids.py, knn2.py); this
    script stays a thin loader.
    """
    usecols = [x_key, y_key, region_key, cluster_col, neighborhood_col]
    dtype = {region_key: "category", cluster_col: "category", neighborhood_col: "category"}
    df = pd.read_csv(raw_csv, usecols=usecols, dtype=dtype)
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

    # Keep the centroid table itself (not just the probabilities scored
    # against it): it is the one number the deleted notebooks stored that is
    # directly, exactly comparable -- see tests/test_parity.py and README's
    # parity table.
    result.uns["neighborhood_centroids"] = centroids

    result.write_h5ad(args.output)


if __name__ == "__main__":
    main()
