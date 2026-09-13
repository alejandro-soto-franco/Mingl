"""Rules for the simulated-transitions parity pipeline (no download needed)."""


rule simulated_transitions_analyse:
    input:
        cells_csv=str(FIXTURES / "synthetic_tissue_{variant}.csv"),
        reference_h5ad=str(FIXTURES / "sim{variant}_results.h5ad"),
    output:
        adata=str(RESULTS / "simulated_transitions" / "{variant}" / "adata.h5ad"),
        parity_json=str(RESULTS / "simulated_transitions" / "{variant}" / "parity.json"),
    params:
        cluster_col=config["simulated_transitions"]["cluster_col"],
        neighborhood_col=config["simulated_transitions"]["neighborhood_col"],
        region_key=config["simulated_transitions"]["region_key"],
        k=config["simulated_transitions"]["k"],
        threshold=config["simulated_transitions"]["threshold"],
    script:
        "../scripts/run_simulated_transitions.py"


rule simulated_transitions_plot:
    input:
        adata=str(RESULTS / "simulated_transitions" / "{variant}" / "adata.h5ad"),
    output:
        pdf=str(RESULTS / "simulated_transitions" / "{variant}" / "figure.pdf"),
        png=str(RESULTS / "simulated_transitions" / "{variant}" / "figure.png"),
    script:
        "../scripts/plot_simulated_transitions.py"
