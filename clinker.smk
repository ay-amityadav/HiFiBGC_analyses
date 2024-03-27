"""
Prepares clinker commands with BGCs in a BiG-SCAPE cluster as input, executes them, and processes the output for further analyses
"""

rule all:
    input:
        expand(f"{{env}}/prepare_clinker_commands.done", env=['sludge', 'chicken', 'human', 'sheep']),
        expand(f"{{env}}/run_clinker.done", env=['sludge', 'chicken', 'human', 'sheep']),
        expand(f"{{env}}/clinker_matrix_lower_mean_analysis.done", env=['sludge', 'chicken', 'human', 'sheep']),

rule prepare_clinker_commands:
    output:
        touch("{env}/prepare_clinker_commands.done")
    run:
        import subprocess
        import pandas as pd

        df_bgc_all_metadata = pd.read_csv(f"HiFiBGC_0.1.13_Run/{wildcards.env}.out/05_final_output/BGC_all_metadata.tsv", sep='\t', index_col=0)

        # create clinker directory
        subprocess.run(["mkdir", "-p", f"{wildcards.env}/clinker"]) 

        clinker_fp = open(f"{wildcards.env}/clinker/clinker_commands.txt", "w")

        # Create clinker command with BGCs belonging to one BiG-SCAPE family 
        df_bgc_all_metadata_grouped = df_bgc_all_metadata.groupby('Family_Number')
        for name, group in df_bgc_all_metadata_grouped:
            name = str(int(name))
            base_dir = f"HiFiBGC_0.1.13_Run/{wildcards.env}.out/04_bgc_clustering/bigscape_input/"
            clinker_command = 'clinker '
            
            subprocess.run(["mkdir", "-p", f"{wildcards.env}/clinker/{name}"])
            
            # get BGC_Id's corresponding to a BiG-SCAPE cluster, here it is index of dataframe group
            for item in group.index:
                file_name = base_dir + item + '.gbk'
                clinker_command += file_name + ' '
            
            # add outputs to clinker_command
            clinker_command += f"-o {wildcards.env}/clinker/{name}/clinker_alignments.csv "
            clinker_command += f"-mo {wildcards.env}/clinker/{name}/matrix_output.csv "
            clinker_command += f"-p {wildcards.env}/clinker/{name}/clinker_plot.html "
            clinker_command += "\n"
            
            # write clinker_command to a file
            clinker_fp.write(clinker_command)

        clinker_fp.close()


rule run_clinker:
    input:
        "{env}/prepare_clinker_commands.done"
    output:
        touch("{env}/run_clinker.done")
    conda:
        "envs/clinker.yml"
    shell:
        """
        while IFS= read -r clinker_command; do

        eval "$clinker_command"

        done < "{wildcards.env}/clinker/clinker_commands.txt"
        """

rule clinker_matrix_lower_mean_analysis:
    input:
        "{env}/run_clinker.done"
    output:
        touch("{env}/clinker_matrix_lower_mean_analysis.done")
    run:
        import glob
        
        import pandas as pd
        import numpy as np

        bigscape_family_number_list = []
        np_clinker_matrix_lower_mean_list = []

        for matrix_file in glob.iglob(f"{wildcards.env}/clinker/*/matrix_output.csv"):
            big_scape_family_number = (matrix_file.split('/'))[-2]
            bigscape_family_number_list.append(big_scape_family_number)
            
            df_clinker_matrix = pd.read_csv(f"{matrix_file}",sep=',', index_col=0)
            np_clinker_matrix = df_clinker_matrix.to_numpy()

            np_clinker_matrix_lower = np_clinker_matrix[np.tril_indices((df_clinker_matrix.shape)[0], -1)]
            np_clinker_matrix_lower_mean = np.mean(np_clinker_matrix_lower)
            
            np_clinker_matrix_lower_mean_list.append(np_clinker_matrix_lower_mean)

        df_mean = pd.DataFrame({'bigscape_family_number': bigscape_family_number_list, 'clinker_lower_matrix_mean': np_clinker_matrix_lower_mean_list})

        # Write to tsv file
        df_mean.to_csv(f"{wildcards.env}/bigscape_family_and_clinker_lower_matrix_mean.tsv", sep='\t', index=False)