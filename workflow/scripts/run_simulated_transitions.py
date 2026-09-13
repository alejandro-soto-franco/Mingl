"""Snakemake script: analyse one SimulatedTransitions variant.

Reproduces the deleted SyntheticGradientTest{Fast,Medium,Slow}.ipynb
notebooks' computation (centroid_Calculation -> cpu_gmm_probability ->
findPositives) on the in-repo synthetic fixture, and reports parity against
the notebooks' own stored output (sim{variant}_results.h5ad).

Snakemake variables (bound by the `script:` directive):
    input.cells_csv, input.reference_h5ad
    output.adata, output.parity_json
    params.cluster_col, params.neighborhood_col, params.region_key,
        params.k, params.threshold
"""

import json

import anndata as ad
import numpy as np
import pandas as pd

from mingl.pp.preprocessing import read_file
from mingl.tl.centroids import centroid_Calculation
from mingl.tl.edges import findPositives
from mingl.tl.gmm import cpu_gmm_probability

snakemake = globals()["snakemake"]

cluster_col = snakemake.params.cluster_col
neighborhood_col = snakemake.params.neighborhood_col
region_key = snakemake.params.region_key
k = int(snakemake.params.k)
threshold = float(snakemake.params.threshold)

cells = read_file(snakemake.input.cells_csv)
cells.obs_names = cells.obs["cell_id"].astype(str)
cells.obs.index.name = None  # avoid an obs index/column name collision on write_h5ad

centroids = centroid_Calculation(
    cells.copy(), k=k, cluster_col=cluster_col, neighborhood_col=neighborhood_col, region_col=region_key
)
result = cpu_gmm_probability(
    cells,
    centroids,
    cluster_col=cluster_col,
    neighborhood_col=neighborhood_col,
    region_key=region_key,
    ks=(k,),
    k=k,
)
findPositives(result, prob_key="neighborhood_probabilities", threshold=threshold)

result.write_h5ad(snakemake.output.adata)

parity = {"variant": snakemake.wildcards.variant, "k": k}
if snakemake.input.reference_h5ad:
    reference = ad.read_h5ad(snakemake.input.reference_h5ad)
    reference_names = list(reference.uns["neighborhood_probability_neighborhoods"])
    reference_probs = pd.DataFrame(
        reference.obsm["neighborhood_probabilities"], index=reference.obs_names, columns=reference_names
    )
    result_names = result.uns["neighborhood_probability_neighborhoods"]
    result_probs = pd.DataFrame(
        result.obsm["neighborhood_probabilities"], index=result.obs_names, columns=result_names
    ).reindex(index=reference_probs.index, columns=reference_probs.columns)

    diff = np.abs(result_probs.values - reference_probs.values)
    parity.update(
        max_abs_diff=float(np.nanmax(diff)),
        mean_abs_diff=float(np.nanmean(diff)),
        pearson_r=float(np.corrcoef(result_probs.values.flatten(), reference_probs.values.flatten())[0, 1]),
        argmax_agreement=float(np.mean(result_probs.idxmax(axis=1).values == reference_probs.idxmax(axis=1).values)),
        reference=str(snakemake.input.reference_h5ad),
    )

with open(snakemake.output.parity_json, "w") as fh:
    json.dump(parity, fh, indent=2)
