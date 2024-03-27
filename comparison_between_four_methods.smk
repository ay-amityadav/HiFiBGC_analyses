"""
This Snakemake file generates dataframes comparing four methods based on total bgc count (further distinguished based on whether partial or complete) and bgc-type count, these dataframes are used for creating visualizations.
"""

rule all:
    input:
        expand("{env}_prepare_bigscape_input.done", env=['sludge', 'chicken', 'human', 'sheep']),
        expand("{env}_run_bigscape.done", env=['sludge', 'chicken', 'human', 'sheep']),
        expand("{env}_four_methods_bgc_bigscape_input_and_output_parse.done", env=['sludge', 'chicken', 'human', 'sheep']),
    
rule prepare_bigscape_input:
    input:
        "HiFiBGC_0.1.13_Run/{env}.out/03_antismash/output/hifiasm-meta",
        "HiFiBGC_0.1.13_Run/{env}.out/03_antismash/output/metaflye",
        "HiFiBGC_0.1.13_Run/{env}.out/03_antismash/output/hicanu",
    output:
        directory("{env}/bigscape_input/hifiasm-meta"),
        directory("{env}/bigscape_input/metaflye"),
        directory("{env}/bigscape_input/hicanu"),
        touch("{env}_prepare_bigscape_input.done"),
    run:
        import glob
        import subprocess
        import os

        for counter in range(len(input)):
            input_dir = input[counter]
            output_dir = output[counter]

            input_dir_files = glob.glob(input_dir + '/*region*.gbk')

            # create output directory
            subprocess.run(["mkdir", "-p", output_dir])

            for file in input_dir_files:
                # get the base name
                file_basename = os.path.basename(file)

                output_file = output_dir + "/" + file_basename

                # create symbolic link
                subprocess.run(["ln", "-rs", file, output_file])


rule install_bigscape:
    output:
        touch("install_bigscape.done"),
    conda:
        "envs/bigscape.yml"
    shell:
        """
        wget https://github.com/medema-group/BiG-SCAPE/archive/refs/tags/v1.1.5.zip
        unzip -o v1.1.5.zip
        rm v1.1.5.zip

        cd BiG-SCAPE-1.1.5
        wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz
        gunzip Pfam-A.hmm.gz
        hmmpress Pfam-A.hmm
        """


rule run_bigscape:
    input:
        "install_bigscape.done",
        bigscape_input_metaflye = "{env}/bigscape_input/metaflye/",
        bigscape_input_hifiasm_meta = "{env}/bigscape_input/hifiasm-meta/",
        bigscape_input_hicanu = "{env}/bigscape_input/hicanu/",
    output:
        touch("{env}_run_bigscape.done"),
        bigscape_out_hifiasm_meta = directory("{env}/bigscape_output/hifiasm-meta"),
        bigscape_out_metaflye = directory("{env}/bigscape_output/metaflye"),
        bigscape_out_hicanu = directory("{env}/bigscape_output/hicanu"),
    conda:
        "envs/bigscape.yml"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {output.bigscape_out_hifiasm_meta}

        python BiG-SCAPE-1.1.5/bigscape.py -i {input.bigscape_input_hifiasm_meta} --cutoffs 0.0001 0.1 0.2 0.3 1.0 \
        --mix --no_classify --hybrids-off --clans-off --cores {threads} -o {output.bigscape_out_hifiasm_meta}

        mkdir -p {output.bigscape_out_metaflye}

        python BiG-SCAPE-1.1.5/bigscape.py -i {input.bigscape_input_metaflye} --cutoffs 0.0001 0.1 0.2 0.3 1.0 \
        --mix --no_classify --hybrids-off --clans-off --cores {threads} -o {output.bigscape_out_metaflye}

        mkdir -p {output.bigscape_out_hicanu}

        python BiG-SCAPE-1.1.5/bigscape.py -i {input.bigscape_input_hicanu} --cutoffs 0.0001 0.1 0.2 0.3 1.0 \
        --mix --no_classify --hybrids-off --clans-off --cores {threads} -o {output.bigscape_out_hicanu}
        """


def four_methods_bgc_bigscape_input_and_output(wildcards):
    import glob
    
    four_methods_bgc_bigscape_input_and_output_list = []

    for assembler in ['metaflye', 'hifiasm-meta', 'hicanu']:
        four_methods_bgc_bigscape_input_and_output_list.append(f"{wildcards.env}/bigscape_input/{assembler}")
    # For hifibgc
    four_methods_bgc_bigscape_input_and_output_list.append(f"HiFiBGC_0.1.13_Run/{wildcards.env}.out/05_final_output/BGC_all")

    for assembler in ['metaflye', 'hifiasm-meta', 'hicanu']:
        for item in glob.iglob(f"{wildcards.env}/bigscape_output/{assembler}/network_files/*/mix/mix_clustering_c0.30.tsv"):
            four_methods_bgc_bigscape_input_and_output_list.append(item)
    # For hifibgc
    for item in glob.iglob(f"HiFiBGC_0.1.13_Run/{wildcards.env}.out/04_bgc_clustering/bigscape_output/network_files/*/mix/mix_clustering_c0.30.tsv"):
        four_methods_bgc_bigscape_input_and_output_list.append(item)

    for assembler in ['metaflye', 'hifiasm-meta', 'hicanu']:
        for item in glob.iglob(f"{wildcards.env}/bigscape_output/{assembler}/network_files/*/Network_Annotations_Full.tsv"):
            four_methods_bgc_bigscape_input_and_output_list.append(item)
    # For hifibgc
    for item in glob.iglob(f"HiFiBGC_0.1.13_Run/{wildcards.env}.out/04_bgc_clustering/bigscape_output/network_files/*/Network_Annotations_Full.tsv"):
        four_methods_bgc_bigscape_input_and_output_list.append(item)

    return(four_methods_bgc_bigscape_input_and_output_list)


rule four_methods_bgc_bigscape_input_and_output_parse:
    input:
        "{env}_run_bigscape.done",
        four_methods_bgc_bigscape_input_and_output,
    output:
        touch("{env}_four_methods_bgc_bigscape_input_and_output_parse.done"),
        "{env}/df_four_methods_representative_bgc_count.tsv",
        "{env}/df_four_methods_bgc_all_metadata.tsv",
    run:
        import os
        import glob
        import subprocess
        import argparse

        from Bio import SeqIO
        import pandas as pd

        import warnings
        warnings.simplefilter(action='ignore', category=FutureWarning)

        def parse_antismash_output(antismash_output_directory):
            """
            Input:
                antismash_output_directory: A folder containing .gbk files that correspond to BGCs predicted by AntiSMASH

            Returns: A dataframe where each row represents a BGC, and the columns hold parsed metadata linked to the respective BGC 
            """
            df_bgcs = pd.DataFrame(columns = ["bgc_id", "contig_id", "length", "contig_edge", "bgc_product"])
            
            for gbk_file in glob.iglob(f"{antismash_output_directory}/*region*.gbk"):
                gbk_file_name = gbk_file.split("/")[-1]
                contig_id = '.'.join((gbk_file_name.split("."))[:-2])
                bgc_id = '.'.join((gbk_file_name.split("."))[:-1])
                records = SeqIO.parse(gbk_file, "gb")

                for record in records:
                    for feature in record.features:
                        qual = feature.qualifiers

                        if feature.type == "region":
                            length = int(feature.location.end) - int(feature.location.start)
                            contig_edge = qual["contig_edge"][0]
                            #bgc_category = qual["category"][0]
                            bgc_product = qual["product"][0]

                            bgc = {"bgc_id": bgc_id, "contig_id": contig_id, "length": length, "contig_edge": contig_edge, "bgc_product": bgc_product}
                            df_bgcs = pd.concat([df_bgcs, pd.DataFrame([bgc])], ignore_index=True)

            return df_bgcs    

        def assign_representative_bgc(df_bgc_merge_bigscape_clustering: pd.DataFrame) -> pd.DataFrame:        
            """
            Assigns a representative BGC among the BGCs in a BiG-SCAPE cluster
            """

            df_bgc_merge_bigscape_clustering_sorted = df_bgc_merge_bigscape_clustering.sort_values(by=['Family Number', 'length', 'contig_edge'], ascending=[True, False, True], na_position='first') # NA's are put at the top of dataframe

            previous_family_number = -1

            for index, row in df_bgc_merge_bigscape_clustering_sorted.iterrows():
                # If 'Family Number' is NA, assign that BGC as representative. 
                if pd.isna(row['Family Number']):
                    df_bgc_merge_bigscape_clustering_sorted.loc[index, 'representative_member'] = True
                    continue

                current_family_number = row['Family Number']
                if (current_family_number != previous_family_number):
                    df_bgc_merge_bigscape_clustering_sorted.loc[index, 'representative_member'] = True
                else:
                    df_bgc_merge_bigscape_clustering_sorted.loc[index, 'representative_member'] = False
                previous_family_number = current_family_number

            return df_bgc_merge_bigscape_clustering_sorted


        df_four_methods_bgc_all_metadata = pd.DataFrame(columns = ['BGC_Id', 'Contig_Id', 'BGC_Length', 'Contig_Edge', 'BGC_Product', 'Family_Number','Product_Prediction','BiG-SCAPE_class', 'Representative_Member', 'Description', 'Organism', 'Taxonomy', 'Method_Name'])
        df_four_methods_representative_bgc_count = pd.DataFrame(columns=['assembler_name', 'representative_bgc_count', 'representative_bgc_partial_count', 'representative_bgc_complete_count'])

        for counter in range(int((len(input[1:]))/3)):        
            bigscape_input_dir = input[1+counter] 
            bigscape_clustering_file = input[1+counter+4] 
            bigscape_network_file = input[1+counter+4+4] 
        
            assembler_name = (bigscape_input_dir.split('/'))[-1]
            if assembler_name == 'BGC_all':
                assembler_name = 'hifibgc'

            # Parse BGC files
            df_bgc = parse_antismash_output(\
                                antismash_output_directory=f"{bigscape_input_dir}")

            # Read BiG-SCAPE output clustering file
            df_bigscape_clustering = pd.read_csv(f"{bigscape_clustering_file}", sep='\t')
            df_bigscape_clustering.columns = ['BGC Name', 'Family Number']

            # Read BiG-SCAPE network annotations file
            df_bigscape_network = pd.read_csv(f"{bigscape_network_file}", sep='\t')
            # columns are: ['BGC', 'Accession ID', 'Description', 'Product Prediction','BiG-SCAPE class', 'Organism', 'Taxonomy'] 

            # In the following statement, 'how='left' is used because BiG-SCAPE occasionally ignores certain BGCs for clustering. Consequently, the size of 'df_bigscape_clustering' may be smaller than that of 'df_bgc'.
            df_bgc_merge_bigscape_clustering = pd.merge(df_bgc, df_bigscape_clustering, how='left', left_on='bgc_id', right_on='BGC Name') 
            
            df_bgc_merge_bigscape_clustering = assign_representative_bgc(df_bgc_merge_bigscape_clustering)
            
            df_bgc_merge_bigscape_clustering_network = pd.merge(df_bgc_merge_bigscape_clustering, df_bigscape_network, how='left', left_on='bgc_id', right_on='BGC') 
            
            # Drop duplicate columns
            df_bgc_merge_bigscape_clustering_network = df_bgc_merge_bigscape_clustering_network.drop(['BGC Name', 'BGC', 'Accession ID'], axis=1)
            # Rename column names
            df_bgc_merge_bigscape_clustering_network.rename(
                columns = {
                    'bgc_id': 'BGC_Id', \
                    'contig_id': 'Contig_Id',
                    'length': 'BGC_Length', 
                    'contig_edge':'Contig_Edge', 
                    'bgc_product': 'BGC_Product',
                    'Family Number': 'Family_Number', 
                    'Product Prediction': 'Product_Prediction', 
                    'BiG-SCAPE class': 'BiG-SCAPE_class',
                    'representative_member': 'Representative_Member'

                }, inplace = True
            ) # there are three other columns: 'Description', 'Organism' and 'Taxonomy'
            
            df_bgc_merge_bigscape_clustering_network['Method_Name'] = assembler_name
            
            df_four_methods_bgc_all_metadata = pd.concat([df_four_methods_bgc_all_metadata, df_bgc_merge_bigscape_clustering_network], ignore_index=True)

            # Subset 'representative_member' BGCs
            df_bgc_merge_bigscape_clustering_network_representative = df_bgc_merge_bigscape_clustering_network[df_bgc_merge_bigscape_clustering_network['Representative_Member'] == True]
            representative_bgc_count = (df_bgc_merge_bigscape_clustering_network_representative.shape)[0]
            representative_bgc_partial_count = (df_bgc_merge_bigscape_clustering_network_representative[df_bgc_merge_bigscape_clustering_network_representative['Contig_Edge'] == "True"]).shape[0]
            representative_bgc_complete_count = (df_bgc_merge_bigscape_clustering_network_representative[df_bgc_merge_bigscape_clustering_network_representative['Contig_Edge'] == "False"]).shape[0]
        
            one_method_representative_bgc_count = pd.DataFrame({'assembler_name': [assembler_name],
                                                                'representative_bgc_count': [representative_bgc_count],
                                                                'representative_bgc_partial_count': [representative_bgc_partial_count],
                                                                'representative_bgc_complete_count': [representative_bgc_complete_count]}
                                                                )
            df_four_methods_representative_bgc_count = pd.concat([df_four_methods_representative_bgc_count, one_method_representative_bgc_count], ignore_index=True)

        df_four_methods_representative_bgc_count.to_csv(output[1], sep='\t', index=False)
        df_four_methods_bgc_all_metadata.to_csv(output[2], sep='\t', index=False)
