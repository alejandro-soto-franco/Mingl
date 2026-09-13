"""Rules for the manuscript's intestine, melanoma and esophagus figures.

Not part of `all` (see Snakefile docstring). Intestine has a confirmed
source (see config.yaml's `datasets.intestine.source`,
`verify_intestine_source` below, and README); melanoma and esophagus do not
and keep their `download_*`/manual-placement path.

Locking: on a machine where the intestine CSV is shared with other
concurrent workflows, set `lock_path` in config (or the MINGL_LOCK_PATH
environment variable) to a lock file and every rule that reads it wraps its
command in `flock` around that path. Neither is set by default, so a fresh
clone or CI runs unlocked, which is correct there: nothing under version
control names a machine-specific path.
"""

import os

_INTESTINE_SOURCE = config["datasets"]["intestine"]["source"]
_LOCK_PATH = config.get("lock_path") or os.environ.get("MINGL_LOCK_PATH")
_LOCK_CMD = f"flock {_LOCK_PATH} " if _LOCK_PATH else ""


def _ks_arg(stage_cfg: dict) -> str:
    return " ".join(str(v) for v in stage_cfg["ks"])


rule verify_intestine_source:
    """Pin the local CSV's identity by SHA-256 before any rule reads it."""
    input:
        raw=_INTESTINE_SOURCE["local_path"],
    output:
        touch(str(RESULTS / "intestine" / ".source_verified")),
    params:
        sha256=_INTESTINE_SOURCE["sha256"],
    shell:
        r"""
        actual=$(sha256sum {input.raw} | cut -d' ' -f1)
        if [ "$actual" != "{params.sha256}" ]; then
            echo "SHA-256 mismatch for {input.raw}: expected {params.sha256}, got $actual" >&2
            exit 1
        fi
        """


rule intestine_neighborhood:
    input:
        raw_csv=_INTESTINE_SOURCE["local_path"],
        verified=str(RESULTS / "intestine" / ".source_verified"),
    output:
        adata=str(RESULTS / "intestine" / "neighborhood" / "adata.h5ad"),
    params:
        cfg=config["datasets"]["intestine"]["neighborhood"],
        ks=_ks_arg(config["datasets"]["intestine"]["neighborhood"]),
    threads: 6
    shell:
        r"""
        nice -n 19 {_LOCK_CMD}python workflow/scripts/run_neighborhood_pipeline.py \
            --raw-csv {input.raw_csv} --output {output.adata} \
            --cluster-col "{params.cfg[cluster_col]}" --neighborhood-col "{params.cfg[neighborhood_col]}" \
            --region-key "{params.cfg[region_key]}" --x-key {params.cfg[x_key]} --y-key {params.cfg[y_key]} \
            --k {params.cfg[k]} --ks {params.ks} --threshold {params.cfg[threshold]}
        """


rule intestine_tissue_unit:
    input:
        raw_csv=_INTESTINE_SOURCE["local_path"],
        verified=str(RESULTS / "intestine" / ".source_verified"),
    output:
        adata=str(RESULTS / "intestine" / "tissue_unit" / "adata.h5ad"),
    params:
        cfg=config["datasets"]["intestine"]["tissue_unit"],
        ks=_ks_arg(config["datasets"]["intestine"]["tissue_unit"]),
    threads: 6
    shell:
        r"""
        nice -n 19 {_LOCK_CMD}python workflow/scripts/run_neighborhood_pipeline.py \
            --raw-csv {input.raw_csv} --output {output.adata} \
            --cluster-col "{params.cfg[cluster_col]}" --neighborhood-col "{params.cfg[neighborhood_col]}" \
            --region-key "{params.cfg[region_key]}" --x-key {params.cfg[x_key]} --y-key {params.cfg[y_key]} \
            --k {params.cfg[k]} --ks {params.ks} --threshold {params.cfg[threshold]}
        """


rule intestine_community:
    input:
        raw_csv=_INTESTINE_SOURCE["local_path"],
        verified=str(RESULTS / "intestine" / ".source_verified"),
    output:
        adata=str(RESULTS / "intestine" / "community" / "adata.h5ad"),
    params:
        cfg=config["datasets"]["intestine"]["community"],
        ks=_ks_arg(config["datasets"]["intestine"]["community"]),
    threads: 6
    shell:
        r"""
        nice -n 19 {_LOCK_CMD}python workflow/scripts/run_neighborhood_pipeline.py \
            --raw-csv {input.raw_csv} --output {output.adata} \
            --cluster-col "{params.cfg[cluster_col]}" --neighborhood-col "{params.cfg[neighborhood_col]}" \
            --region-key "{params.cfg[region_key]}" --x-key {params.cfg[x_key]} --y-key {params.cfg[y_key]} \
            --k {params.cfg[k]} --ks {params.ks} --threshold {params.cfg[threshold]}
        """


rule intestine_networks:
    """fig3: interaction-pair graphs at all three hierarchical levels.

    Reuses the three adata.h5ad above (already-computed probabilities), so
    it reads no large CSV and needs no lock.
    """
    input:
        neighborhood=str(RESULTS / "intestine" / "neighborhood" / "adata.h5ad"),
        tissue_unit=str(RESULTS / "intestine" / "tissue_unit" / "adata.h5ad"),
        community=str(RESULTS / "intestine" / "community" / "adata.h5ad"),
    output:
        summary=str(RESULTS / "intestine" / "networks" / "summary.json"),
        neighborhood_pdf=str(RESULTS / "intestine" / "networks" / "neighborhood.pdf"),
        neighborhood_png=str(RESULTS / "intestine" / "networks" / "neighborhood.png"),
        tissue_unit_pdf=str(RESULTS / "intestine" / "networks" / "tissue_unit.pdf"),
        tissue_unit_png=str(RESULTS / "intestine" / "networks" / "tissue_unit.png"),
        community_pdf=str(RESULTS / "intestine" / "networks" / "community.pdf"),
        community_png=str(RESULTS / "intestine" / "networks" / "community.png"),
    params:
        threshold=config["datasets"]["intestine"]["neighborhood"]["threshold"],
    script:
        "../scripts/run_intestine_networks.py"


rule download_melanoma:
    output:
        raw=str(DATA / "melanoma" / "raw.csv"),
    params:
        url=config["datasets"]["melanoma"]["source"]["download_url"],
        sha256=config["datasets"]["melanoma"]["source"]["sha256"],
        description=config["datasets"]["melanoma"]["source"]["description"],
    script:
        "../scripts/fetch_dataset.py"


rule melanoma_neighborhood:
    input:
        raw_csv=str(DATA / "melanoma" / "raw.csv"),
    output:
        adata=str(RESULTS / "melanoma" / "neighborhood" / "adata.h5ad"),
    params:
        cfg=config["datasets"]["melanoma"]["neighborhood"],
        ks=_ks_arg(config["datasets"]["melanoma"]["neighborhood"]),
    threads: 6
    shell:
        r"""
        nice -n 19 {_LOCK_CMD}python workflow/scripts/run_neighborhood_pipeline.py \
            --raw-csv {input.raw_csv} --output {output.adata} \
            --cluster-col "{params.cfg[cluster_col]}" --neighborhood-col "{params.cfg[neighborhood_col]}" \
            --region-key "{params.cfg[region_key]}" --x-key {params.cfg[x_key]} --y-key {params.cfg[y_key]} \
            --k {params.cfg[k]} --ks {params.ks} --threshold {params.cfg[threshold]}
        """


rule esophagus_neighborhood:
    input:
        raw_csv=str(DATA / "esophagus" / "raw.csv"),
    output:
        adata=str(RESULTS / "esophagus" / "neighborhood" / "adata.h5ad"),
    params:
        cfg=config["datasets"]["esophagus"]["neighborhood"],
        ks=_ks_arg(config["datasets"]["esophagus"]["neighborhood"]),
    threads: 6
    shell:
        r"""
        nice -n 19 python workflow/scripts/run_neighborhood_pipeline.py \
            --raw-csv {input.raw_csv} --output {output.adata} \
            --cluster-col "{params.cfg[cluster_col]}" --neighborhood-col "{params.cfg[neighborhood_col]}" \
            --region-key "{params.cfg[region_key]}" --k {params.cfg[k]} --ks {params.ks} \
            --threshold {params.cfg[threshold]}
        """
