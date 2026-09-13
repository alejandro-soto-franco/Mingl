"""Regression test for crd()'s configurable id/region/neighbourhood columns.

Previously ``crd`` hardcoded the literals 'cellid', 'region' and 'neigh_name'
throughout its body; this checks the same computation still runs, and now
also runs, under caller-supplied column names.
"""

import anndata as ad
import numpy as np
import pandas as pd

from mingl.tl.crd import crd


def _small_case(cellid_col: str, region_col: str, neigh_name_col: str):
    obs = pd.DataFrame(
        {
            neigh_name_col: ["N1", "N1", "N2", "N2"],
            region_col: ["R1", "R1", "R1", "R1"],
        },
        index=["c1", "c2", "c3", "c4"],
    )
    obs[cellid_col] = obs.index
    cells = ad.AnnData(X=np.zeros((4, 0)), obs=obs)

    windows2 = pd.DataFrame(
        {
            cellid_col: ["c1", "c2", "c3", "c4"],
            region_col: ["R1", "R1", "R1", "R1"],
            "A": [1.0, 0.8, 0.1, 0.2],
            "B": [0.0, 0.2, 0.9, 0.8],
        },
        index=["c1", "c2", "c3", "c4"],
    )

    probabilities_df = pd.DataFrame(
        {
            cellid_col: ["c1", "c2", "c3", "c4"],
            "N1": [0.6, 0.6, 0.4, 0.4],
            "N2": [0.4, 0.4, 0.6, 0.6],
        }
    )
    return cells, windows2, probabilities_df


def test_crd_runs_with_default_column_names(tmp_path):
    cells, windows2, probabilities_df = _small_case("cellid", "region", "neigh_name")

    probs_df, deltas_df = crd(
        cells,
        windows2,
        probabilities_df,
        cell_type_features=["A", "B"],
        out_probs_path=str(tmp_path / "probs.csv"),
        out_delta_path=str(tmp_path / "deltas.csv"),
    )

    assert {"N1", "N2", "cellid", "region", "neigh_name"}.issubset(probs_df.columns)
    assert len(probs_df) == 4


def test_crd_runs_with_custom_column_names(tmp_path):
    cells, windows2, probabilities_df = _small_case("cell_id", "tissue_region", "assigned_neighborhood")

    probs_df, deltas_df = crd(
        cells,
        windows2,
        probabilities_df,
        cell_type_features=["A", "B"],
        cellid_col="cell_id",
        region_col="tissue_region",
        neigh_name_col="assigned_neighborhood",
        out_probs_path=str(tmp_path / "probs.csv"),
        out_delta_path=str(tmp_path / "deltas.csv"),
    )

    assert {"N1", "N2", "cell_id", "tissue_region", "assigned_neighborhood"}.issubset(probs_df.columns)
    assert len(probs_df) == 4
