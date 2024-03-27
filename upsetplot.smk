"""
Snakemake for creating upsetplots
"""

import glob

C_PARAMETER = "c0.30"

rule all:
    input:  
        expand(f"{{env}}/df_upsetplot.tsv", env=['sludge', 'human', 'sheep', 'chicken']),
        expand(f"{{env}}/upsetplot_{C_PARAMETER}.pdf", env=['sludge', 'human', 'sheep', 'chicken']),
        expand(f"{{env}}/upsetplot_{C_PARAMETER}_with_percentages.pdf", env=['sludge', 'human', 'sheep', 'chicken']),
        expand(f"{{env}}/upsetplot_{C_PARAMETER}_clinker_lower_matrix_mean.pdf", env=['sludge', 'human', 'sheep', 'chicken']),
        expand(f"{{env}}/upsetplot_{C_PARAMETER}_bigscape_class.pdf", env=['sludge', 'human', 'sheep', 'chicken']),

rule upsetplot_prepare_input:
    input:
        "HiFiBGC_0.1.13_Run/{env}.out/05_final_output/upsetplot/df_bgc_family_to_dataset_c0.30.tsv",
        "HiFiBGC_0.1.13_Run/{env}.out/05_final_output/BGC_all_metadata.tsv",
        "{env}/bigscape_family_and_clinker_lower_matrix_mean.tsv",
    output:
        "{env}/df_upsetplot.tsv",
    shell:
        """
        mkdir -p {wildcards.env}

        python3 upsetplot_prepare_input.py \
	        --df_bgc_family_to_dataset {input[0]} \
	        --bgc_all_metadata  {input[1]} \
            --bigscape_family_and_clinker_lower_matrix_mean {input[2]} \
	        --df_upsetplot {output}
        """

rule plot_upsetplot:
    input:
        "{env}/df_upsetplot.tsv",
    output:
        f"{{env}}/upsetplot_{C_PARAMETER}.pdf",
        f"{{env}}/upsetplot_{C_PARAMETER}_with_percentages.pdf",
        f"{{env}}/upsetplot_{C_PARAMETER}_clinker_lower_matrix_mean.pdf",
        f"{{env}}/upsetplot_{C_PARAMETER}_bigscape_class.pdf",
    conda:
        "envs/r_complexupset.yml"
    shell:
        """
        mkdir -p {wildcards.env}

        Rscript upsetplot.R {input} {wildcards.env} {C_PARAMETER}
        """

