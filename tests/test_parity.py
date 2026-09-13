"""End-to-end parity: re-run the neighbourhood-probability pipeline on the
shipped SimulatedTransitions fixtures and compare against the notebooks'
stored ``sim*_results.h5ad`` outputs.

These three h5ad files are exactly what the deleted
``SyntheticGradientTest{Fast,Medium,Slow}.ipynb`` notebooks produced (seed=42,
k=10 throughout, per ``FunctionsSimulator.ipynb``): each stores the per-cell
input table plus the centroid AnnData (``uns["neighborhood_centroids"]``) and
neighbourhood-probability matrix (``obsm["neighborhood_probabilities"]``) the
notebooks computed via ``centroid_Calculation`` + ``cpu_gmm_probability``.

Feeding the shipped centroids back through ``cpu_gmm_probability`` at k=10
reproduces the stored probabilities exactly, which is this repository's
strongest available parity evidence: unlike the manuscript's intestine,
melanoma and esophagus figures (whose source data this workflow could not
obtain in a schema-confirmed form; see README's "Differences from upstream"),
this fixture is a real notebook output already in the repository, with no
external dependency and no ambiguity about what pipeline produced it.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from mingl.pp.preprocessing import read_file
from mingl.tl.centroids import centroid_Calculation
from mingl.tl.gmm import cpu_gmm_probability

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "simulated_transitions"


@pytest.mark.parametrize("name", ["fast", "medium", "slow"])
def test_neighborhood_probabilities_match_shipped_notebook_output(name):
    reference = ad.read_h5ad(FIXTURE_DIR / f"sim{name}_results.h5ad")
    reference_names = list(reference.uns["neighborhood_probability_neighborhoods"])
    reference_probs = pd.DataFrame(
        reference.obsm["neighborhood_probabilities"], index=reference.obs_names, columns=reference_names
    )

    cells = read_file(FIXTURE_DIR / f"synthetic_tissue_{name}.csv")
    cells.obs_names = cells.obs["cell_id"].astype(str)

    centroids = centroid_Calculation(
        cells.copy(), k=10, cluster_col="cell_type", neighborhood_col="Neighborhood", region_col="unique_region"
    )
    result = cpu_gmm_probability(
        cells.copy(),
        centroids,
        cluster_col="cell_type",
        neighborhood_col="Neighborhood",
        region_key="unique_region",
        ks=(10,),
        k=10,
    )
    result_names = result.uns["neighborhood_probability_neighborhoods"]
    result_probs = pd.DataFrame(
        result.obsm["neighborhood_probabilities"], index=result.obs_names, columns=result_names
    ).reindex(index=reference_probs.index, columns=reference_probs.columns)

    np.testing.assert_allclose(result_probs.values, reference_probs.values, atol=1e-9)
