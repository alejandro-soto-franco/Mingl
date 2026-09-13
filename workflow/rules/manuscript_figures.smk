"""Rules for the manuscript's intestine, melanoma and esophagus figures.

Not part of `all` (see Snakefile docstring): these need raw data this
repository does not ship. Each dataset's `download_*` rule fetches and
SHA-256-verifies its documented source (config.yaml `datasets.<name>.source`)
into DATA/<name>/raw.csv; esophagus has no download rule because no public
source was identified (see README) -- place a file at
`data/esophagus/raw.csv` yourself to run it.
"""


rule download_intestine:
    output:
        raw=str(DATA / "intestine" / "raw.csv"),
    params:
        url=config["datasets"]["intestine"]["source"]["download_url"],
        sha256=config["datasets"]["intestine"]["source"]["sha256"],
        description=config["datasets"]["intestine"]["source"]["description"],
    script:
        "../scripts/fetch_dataset.py"


rule download_melanoma:
    output:
        raw=str(DATA / "melanoma" / "raw.csv"),
    params:
        url=config["datasets"]["melanoma"]["source"]["download_url"],
        sha256=config["datasets"]["melanoma"]["source"]["sha256"],
        description=config["datasets"]["melanoma"]["source"]["description"],
    script:
        "../scripts/fetch_dataset.py"


rule intestine_neighborhood:
    input:
        raw_csv=str(DATA / "intestine" / "raw.csv"),
    output:
        adata=str(RESULTS / "intestine" / "neighborhood" / "adata.h5ad"),
    params:
        **config["datasets"]["intestine"]["neighborhood"],
    script:
        "../scripts/run_neighborhood_pipeline.py"


rule intestine_tissue_unit:
    input:
        raw_csv=str(DATA / "intestine" / "raw.csv"),
    output:
        adata=str(RESULTS / "intestine" / "tissue_unit" / "adata.h5ad"),
    params:
        **config["datasets"]["intestine"]["tissue_unit"],
    script:
        "../scripts/run_neighborhood_pipeline.py"


rule intestine_community:
    input:
        raw_csv=str(DATA / "intestine" / "raw.csv"),
    output:
        adata=str(RESULTS / "intestine" / "community" / "adata.h5ad"),
    params:
        **config["datasets"]["intestine"]["community"],
    script:
        "../scripts/run_neighborhood_pipeline.py"


rule melanoma_neighborhood:
    input:
        raw_csv=str(DATA / "melanoma" / "raw.csv"),
    output:
        adata=str(RESULTS / "melanoma" / "neighborhood" / "adata.h5ad"),
    params:
        **config["datasets"]["melanoma"]["neighborhood"],
    script:
        "../scripts/run_neighborhood_pipeline.py"


rule esophagus_neighborhood:
    input:
        raw_csv=str(DATA / "esophagus" / "raw.csv"),
    output:
        adata=str(RESULTS / "esophagus" / "neighborhood" / "adata.h5ad"),
    params:
        cluster_col=config["datasets"]["esophagus"]["neighborhood"]["cluster_col"],
        neighborhood_col=config["datasets"]["esophagus"]["neighborhood"]["neighborhood_col"],
        region_key=config["datasets"]["esophagus"]["neighborhood"]["region_key"],
        k=config["datasets"]["esophagus"]["neighborhood"]["k"],
        ks=config["datasets"]["esophagus"]["neighborhood"]["ks"],
        threshold=config["datasets"]["esophagus"]["neighborhood"]["threshold"],
    script:
        "../scripts/run_neighborhood_pipeline.py"
