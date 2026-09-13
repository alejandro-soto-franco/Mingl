"""Snakemake script: generic k-window neighbourhood/tissue-unit/community pass.

The manuscript's intestine (fig2, fig4) and melanoma (fig2) figures share one
computation once their raw table is loaded: KNN2 windows -> per-neighbourhood
centroids -> Gaussian-mixture membership probabilities -> threshold count.
`centroid_Calculation` already calls KNN2 internally, so this script only
calls it and `cpu_gmm_probability`, parameterized entirely from config.yaml
(`datasets.<name>.neighborhood` / `.tissue_unit` / `.community`) rather than
the manuscript-specific defaults the library used to hardcode.

Snakemake variables:
    input.raw_csv
    output.adata
    params.cluster_col, params.neighborhood_col, params.region_key,
        params.k, params.ks, params.threshold
"""

from mingl.pp.preprocessing import read_file
from mingl.tl.centroids import centroid_Calculation
from mingl.tl.edges import findPositives
from mingl.tl.gmm import cpu_gmm_probability

snakemake = globals()["snakemake"]

cluster_col = snakemake.params.cluster_col
neighborhood_col = snakemake.params.neighborhood_col
region_key = snakemake.params.region_key
k = int(snakemake.params.k)
ks = tuple(int(v) for v in snakemake.params.ks)
threshold = float(snakemake.params.threshold)

cells = read_file(snakemake.input.raw_csv)

centroids = centroid_Calculation(
    cells.copy(), k=k, cluster_col=cluster_col, neighborhood_col=neighborhood_col, region_col=region_key
)
result = cpu_gmm_probability(
    cells,
    centroids,
    cluster_col=cluster_col,
    neighborhood_col=neighborhood_col,
    region_key=region_key,
    ks=ks,
    k=k,
)
findPositives(result, prob_key="neighborhood_probabilities", threshold=threshold)

result.write_h5ad(snakemake.output.adata)
