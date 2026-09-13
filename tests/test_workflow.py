"""Dry-run checks that the Snakemake workflow's DAG is well-formed.

These do not execute any rule (``-n``/``--dry-run``); they only confirm
Snakemake can parse workflow/Snakefile against both configs and resolve a
concrete plan for the default target. Actually running the pipeline is
covered by ``pixi run smoke`` (config/smoke.yaml, real computation on the
in-repo synthetic fixture) and by tests/test_parity.py.
"""

import shutil
import subprocess
import sys

import pytest

pytestmark = pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not on PATH")


def _dry_run(*extra_args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "--snakefile",
            "workflow/Snakefile",
            "--configfile",
            "config/config.yaml",
            *extra_args,
            "-n",
            "--forceall",
        ],
        capture_output=True,
        text=True,
    )


def test_dry_run_default_config():
    result = _dry_run()
    assert result.returncode == 0, result.stderr


def test_dry_run_smoke_config():
    result = _dry_run("config/smoke.yaml")
    assert result.returncode == 0, result.stderr
    # the smoke config restricts simulated_transitions to one variant
    assert "fast" in result.stdout
