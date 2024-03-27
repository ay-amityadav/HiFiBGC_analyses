"""
Prepares input file for r-complexupset package for creating upset plots
"""

import argparse
import glob

import pandas as pd

def main():
    # Argument parsing
    parser = argparse.ArgumentParser()

    parser.add_argument("--df_bgc_family_to_dataset", type=str, required=True)
    parser.add_argument("--bgc_all_metadata", type=str, required=True)
    parser.add_argument("--bigscape_family_and_clinker_lower_matrix_mean", type=str, required=True)
    parser.add_argument("--df_upsetplot", type=str, required=True)

    args = parser.parse_args()

    df_bgc_family_to_dataset = args.df_bgc_family_to_dataset
    bgc_all_metadata = args.bgc_all_metadata
    bigscape_family_and_clinker_lower_matrix_mean = args.bigscape_family_and_clinker_lower_matrix_mean
    df_upsetplot = args.df_upsetplot

    df_bigscape_upset = pd.read_csv(f"{df_bgc_family_to_dataset}", sep='\t', index_col=0, dtype={'hicanu': str, 'metaflye': str, 'hifiasm-meta': str, 'unmapped_reads': str, 'Family_Number': str})

    df_bgc_metadata = pd.read_csv(f"{bgc_all_metadata}", sep='\t')
    df_bgc_metadata_representative = df_bgc_metadata[df_bgc_metadata['Representative_Member'] == True]
    df_bgc_metadata_representative = df_bgc_metadata_representative[['BGC_Id', 'Family_Number', 'Contig_Edge', 'Representative_Member', 'BGC_Product']]

    # For BGCs with Family_Number=NA, set Family_Number to BGC_Id; this is done for merging dataframes in below section
    for (index, row) in df_bgc_metadata_representative.iterrows():
        if not pd.isna(row['Family_Number']):
            df_bgc_metadata_representative.loc[index, 'Family_Number'] = str(int(df_bgc_metadata_representative.loc[index, 'Family_Number']))
        else:
            df_bgc_metadata_representative.loc[index, 'Family_Number'] = df_bgc_metadata_representative.loc[index, 'BGC_Id']

    # Add BGC metadata to upsetplot
    assert (df_bigscape_upset.shape)[0] == (df_bgc_metadata_representative.shape)[0]
    df_bigscape_upset = pd.merge(df_bigscape_upset, df_bgc_metadata_representative, on='Family_Number', how='inner')
    assert (df_bigscape_upset.shape)[0] == (df_bgc_metadata_representative.shape)[0]


    # Add clinker_lower_matrix_mean component to upsetplot
    df_bigscape_family_and_clinker_lower_matrix_mean = pd.read_csv(f"{bigscape_family_and_clinker_lower_matrix_mean}", sep='\t')
    df_bigscape_family_and_clinker_lower_matrix_mean['bigscape_family_number'] = df_bigscape_family_and_clinker_lower_matrix_mean['bigscape_family_number'].apply(lambda x: str(x))
    df_bigscape_upset = pd.merge(df_bigscape_upset, df_bigscape_family_and_clinker_lower_matrix_mean, how='left', left_on='Family_Number', right_on='bigscape_family_number') # how='outer' can also be used.


    # Add `BiG-SCAPE class` component to upsetplot
    env = ((bgc_all_metadata.split('/')[1]).split('.'))[0]
    x = glob.iglob(f"HiFiBGC_0.1.13_Run/{env}.out/04_bgc_clustering/bigscape_output/network_files/*/Network_Annotations_Full.tsv")
    bigscape_network_annotations_file = (list(x))[0] # there is only one item in list
    df_bigscape_network_annotations = pd.read_csv(f"{bigscape_network_annotations_file}", sep='\t')
    df_bigscape_upset = pd.merge(df_bigscape_upset, df_bigscape_network_annotations, how='left', left_on='BGC_Id', right_on='BGC')


    # For BGCs ignored by BiG-SCAPE, 'BiG-SCAPE class' column is returned as NA by BiG-SCAPE
    # On inspecting other columns associated with such BGCs, these were found to be `RiPPs`, across datasets having any such BGCs
    df_bigscape_upset['BiG-SCAPE class'] = df_bigscape_upset['BiG-SCAPE class'].fillna('RiPPs')

    # Select columns
    if 'unmapped_reads' in df_bigscape_upset.columns:
        df_bigscape_upset = df_bigscape_upset[["Family_Number", "hicanu", "metaflye", "hifiasm-meta", "unmapped_reads", "Contig_Edge", "clinker_lower_matrix_mean", "BiG-SCAPE class"]]    
    else:
        df_bigscape_upset = df_bigscape_upset[["Family_Number", "hicanu", "metaflye", "hifiasm-meta", "Contig_Edge", "clinker_lower_matrix_mean", "BiG-SCAPE class"]]

    # Rename column name
    df_bigscape_upset.rename(columns = {'Contig_Edge':'Partial/Complete', 'BiG-SCAPE class':'BiGSCAPE_class'}, inplace = True) 

    # Write out to a file
    df_bigscape_upset.to_csv(f"{df_upsetplot}", sep='\t', index=False)

if __name__ == '__main__':
    main()