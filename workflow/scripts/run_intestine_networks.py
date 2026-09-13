"""Snakemake script: fig3, interaction-pair graphs at three hierarchical levels.

Reads the three already-computed level adata.h5ad files (neighbourhood,
tissue-unit, community; each has its neighbourhood-probability columns
copied into .obs by run_neighborhood_pipeline.py) and, for each, builds and
plots the top-15 neighbourhood-pair interaction graph.

Snakemake variables:
    input.neighborhood, input.tissue_unit, input.community
    output.summary, output.{neighborhood,tissue_unit,community}_{pdf,png}
    params.threshold
"""

import json

import anndata as ad
import matplotlib

matplotlib.use("Agg")
# Node labels here are raw neighbourhood/community/tissue-unit names from the
# data (e.g. "Stroma & Innate Immune"), not typeset mathematical notation;
# LaTeX treats '&' as a table column separator and fails outright on it.
# Scoped to this script only, not a global rcParams/matplotlibrc change.
matplotlib.rcParams["text.usetex"] = False

from mingl.tl.network_graphs import build_neighborhood_pair_graph, plot_neighborhood_pair_graph  # noqa: E402

snakemake = globals()["snakemake"]

threshold = float(snakemake.params.threshold)
levels = {
    "neighborhood": (
        snakemake.input.neighborhood,
        snakemake.output.neighborhood_pdf,
        snakemake.output.neighborhood_png,
    ),
    "tissue_unit": (snakemake.input.tissue_unit, snakemake.output.tissue_unit_pdf, snakemake.output.tissue_unit_png),
    "community": (snakemake.input.community, snakemake.output.community_pdf, snakemake.output.community_png),
}

summary = {}
for level, (adata_path, pdf_path, png_path) in levels.items():
    adata = ad.read_h5ad(adata_path)
    prob_cols = list(adata.uns["neighborhood_probability_neighborhoods"])

    graph, top_pairs = build_neighborhood_pair_graph(adata, prob_cols, threshold=threshold)
    fig = plot_neighborhood_pair_graph(adata, return_fig=True)
    fig.savefig(pdf_path)
    fig.savefig(png_path)

    summary[level] = {
        "n_neighborhoods": len(prob_cols),
        "n_edges": graph.number_of_edges(),
        "top_pairs": top_pairs.to_dict(orient="records"),
    }

with open(snakemake.output.summary, "w") as fh:
    json.dump(summary, fh, indent=2, default=str)
