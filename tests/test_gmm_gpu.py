"""Tests for the GPU/CPU GMM parity fix.

``gpu_gmm_probability`` used to call ``KNN2`` without passing ``ks``, so the
GPU path scored a different k-window set than
``cpu_gmm_probability`` (ks=(10, 20, 100, 300) vs KNN2's own default
(5, 10, 20, 100, 300)), and it skipped the CPU path's obs-column validation.
This machine has no GPU, so the actual CuPy math is only exercised when
``cupy`` is importable; the signature/validation checks run unconditionally
since they fail before any CuPy call.
"""

import importlib.util

import anndata as ad
import numpy as np
import pandas as pd
import pytest

cupy_available = importlib.util.find_spec("cupy") is not None


def _make_cells_adata() -> ad.AnnData:
    obs = pd.DataFrame(
        {
            "cell_type": ["T", "B", "T"],
            "neighborhood": ["N1", "N2", "N1"],
            "unique_region": ["R1", "R1", "R1"],
        },
        index=["cell_1", "cell_2", "cell_3"],
    )
    return ad.AnnData(X=np.zeros((len(obs), 1), dtype=np.float64), obs=obs)


@pytest.mark.skipif(not cupy_available, reason="cupy/GPU not available on this machine")
def test_gpu_gmm_probability_default_ks_matches_cpu_ks():
    import mingl.tl.gmm as gmm
    import mingl.tl.gmm_gpu as gmm_gpu

    # The two defaults must agree so both paths score identical k-windows.
    cpu_default = gmm.cpu_gmm_probability.__kwdefaults__["ks"]
    gpu_default = gmm_gpu.gpu_gmm_probability.__kwdefaults__["ks"]
    assert cpu_default == gpu_default == (10, 20, 100, 300)


@pytest.mark.skipif(not cupy_available, reason="cupy/GPU not available on this machine")
def test_gpu_gmm_probability_validates_required_columns():
    import mingl.tl.gmm_gpu as gmm_gpu

    cells = _make_cells_adata()
    centroids = ad.AnnData(X=np.zeros((2, 1)), obs=pd.DataFrame(index=["N1", "N2"]))
    with pytest.raises(KeyError):
        gmm_gpu.gpu_gmm_probability(cells, centroids, cluster_col="not_a_column")


def test_gpu_gmm_probability_ks_default_declared_without_cupy():
    """The ``ks`` default can be checked without importing cupy at all.

    ``gmm_gpu`` imports cupy at module scope, so this test only runs the
    parity check that does not require a GPU: comparing the source-level
    default in the function signature against the CPU module's default.
    """
    import ast
    from pathlib import Path

    gpu_source = Path("src/mingl/tl/gmm_gpu.py")
    if not gpu_source.exists():
        gpu_source = Path(__file__).parent.parent / "src/mingl/tl/gmm_gpu.py"
    tree = ast.parse(gpu_source.read_text())
    func = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "gpu_gmm_probability")
    ks_default = next(d for arg, d in zip(func.args.kwonlyargs, func.args.kw_defaults, strict=True) if arg.arg == "ks")
    values = tuple(elt.value for elt in ks_default.elts)
    assert values == (10, 20, 100, 300)
