"""Snakemake script: plot one SimulatedTransitions variant's neighbourhood map.

Reads the AnnData written by run_simulated_transitions.py and draws a
spatial scatter of cells coloured by their most probable neighbourhood,
matching the "figures are rule outputs (PDF + PNG)" convention.

Snakemake variables:
    input.adata
    output.pdf, output.png
    wildcards.variant
"""

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

snakemake = globals()["snakemake"]

adata = ad.read_h5ad(snakemake.input.adata)
names = adata.uns["neighborhood_probability_neighborhoods"]
probs = pd.DataFrame(adata.obsm["neighborhood_probabilities"], index=adata.obs_names, columns=names)
dominant = probs.idxmax(axis=1)

fig, ax = plt.subplots(figsize=(6, 6), dpi=200)
for name, group in dominant.groupby(dominant):
    idx = group.index
    ax.scatter(adata.obs.loc[idx, "x"], adata.obs.loc[idx, "y"], s=6, label=str(name), alpha=0.8)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title(f"Simulated transitions ({snakemake.wildcards.variant})")
ax.set_aspect("equal", "box")
ax.legend(title="Neighborhood", frameon=False, markerscale=2)
fig.tight_layout()

fig.savefig(snakemake.output.pdf)
fig.savefig(snakemake.output.png)
plt.close(fig)
