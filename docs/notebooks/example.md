---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  name: python3
---

# Example

A minimal, runnable walk through the neighbourhood-probability pipeline this
package implements: k-window neighbourhood composition, per-neighbourhood
centroids, and Gaussian-mixture membership probabilities. It uses a tiny
synthetic dataset built in place, so it executes on every documentation
build with no external data.

```{code-cell} ipython3
import numpy as np
import pandas as pd
import anndata as ad

import mingl
```

Build a synthetic set of cells: two regions, two cell types, one
pre-assigned neighbourhood label per cell.

```{code-cell} ipython3
rng = np.random.default_rng(0)
# KNN2 defaults to k-windows up to 300 neighbours per region, so each region
# needs more than 300 cells.
n = 800
obs = pd.DataFrame(
    {
        "x": rng.uniform(0, 100, n),
        "y": rng.uniform(0, 100, n),
        "unique_region": rng.choice(["R1", "R2"], n),
        "cell_type": rng.choice(["A", "B"], n),
        "neighborhood": rng.choice(["N1", "N2"], n),
    }
)
adata = ad.AnnData(X=np.zeros((n, 0), dtype=np.float32), obs=obs)
adata
```

Compute per-neighbourhood centroids from k=10 windows, then score every
cell's membership probability against each neighbourhood.

```{code-cell} ipython3
centroids = mingl.tl.centroid_Calculation(
    adata, k=10, cluster_col="cell_type", neighborhood_col="neighborhood", region_col="unique_region"
)
result = mingl.tl.cpu_gmm_probability(
    adata, centroids, cluster_col="cell_type", neighborhood_col="neighborhood",
    region_key="unique_region", ks=(10,), k=10,
)
pd.DataFrame(
    result.obsm["neighborhood_probabilities"],
    columns=result.uns["neighborhood_probability_neighborhoods"],
).head()
```

See {func}`mingl.tl.centroid_Calculation`, {func}`mingl.tl.gmm.cpu_gmm_probability`
and {func}`mingl.tl.KNN2` for the full parameter reference, and
`workflow/Snakefile` for how this same pipeline reproduces the manuscript
figures end to end.
