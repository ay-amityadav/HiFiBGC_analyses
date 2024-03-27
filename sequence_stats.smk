
"""
This snakemake computes:
    1. Fraction of count of unmapped reads to count of initial raw reads
    2. N50 statistic for the three assemblies
"""

rule all:
    input:
        expand(f"{{env}}/sequence_stats.done", env=['sludge', 'chicken', 'human', 'sheep']),

rule seqkit_stats:
    input:  
        raw_reads="raw/{env}.fastq",
        unmapped_reads="HiFiBGC_0.1.13_Run/{env}.out/02_mapping_reads_to_merged_assembly/reads_mapped_to_merged_assembly_unmapped.fasta",
        hifiasm_meta_assembly="HiFiBGC_0.1.13_Run/{env}.out/01_assembly/hifiasm-meta/hifiasm_meta.p_contigs.fa",
        metaflye_assembly="HiFiBGC_0.1.13_Run/{env}.out/01_assembly/metaflye/assembly.fasta",
        hicanu_assembly="HiFiBGC_0.1.13_Run/{env}.out/01_assembly/hicanu/hicanu.contigs.fasta",
    output:
        raw_reads_seqkit_stats=temp("{env}/tmp_raw_reads_seqkit_stats.tsv"),
        unmapped_reads_seqkit_stats=temp("{env}/tmp_unmapped_reads_seqkit_stats.tsv"),
        hifiasm_meta_assembly_seqkit_stats=temp("{env}/tmp_hifiasm_meta_assembly_seqkit_stats.tsv"),
        metaflye_assembly_seqkit_stats=temp("{env}/tmp_metaflye_assembly_seqkit_stats.tsv"),
        hicanu_assembly_seqkit_stats=temp("{env}/tmp_hicanu_assembly_seqkit_stats.tsv"),
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit stats -a -T {input.raw_reads} -o {output.raw_reads_seqkit_stats}

        seqkit stats -a -T {input.unmapped_reads} -o {output.unmapped_reads_seqkit_stats}

        seqkit stats -a -T {input.hifiasm_meta_assembly} -o {output.hifiasm_meta_assembly_seqkit_stats}

        seqkit stats -a -T {input.metaflye_assembly} -o {output.metaflye_assembly_seqkit_stats}

        seqkit stats -a -T {input.hicanu_assembly} -o {output.hicanu_assembly_seqkit_stats}
        """

rule get_statistics:
    input:
        raw_reads_seqkit_stats="{env}/tmp_raw_reads_seqkit_stats.tsv",
        unmapped_reads_seqkit_stats="{env}/tmp_unmapped_reads_seqkit_stats.tsv",
        hifiasm_meta_assembly_seqkit_stats="{env}/tmp_hifiasm_meta_assembly_seqkit_stats.tsv",
        metaflye_assembly_seqkit_stats="{env}/tmp_metaflye_assembly_seqkit_stats.tsv",
        hicanu_assembly_seqkit_stats="{env}/tmp_hicanu_assembly_seqkit_stats.tsv",
    output:
        touch("{env}/sequence_stats.done"),
    run:
       import pandas as pd

       df_raw_reads_seqkit_stats = pd.read_csv(f"{input.raw_reads_seqkit_stats}", sep='\t')
       df_unmapped_reads_seqkit_stats = pd.read_csv(f"{input.unmapped_reads_seqkit_stats}", sep='\t')
       df_hifiasm_meta_assembly_seqkit_stats = pd.read_csv(f"{input.hifiasm_meta_assembly_seqkit_stats}", sep='\t')
       df_metaflye_assembly_seqkit_stats = pd.read_csv(f"{input.metaflye_assembly_seqkit_stats}", sep='\t')
       df_hicanu_assembly_seqkit_stats = pd.read_csv(f"{input.hicanu_assembly_seqkit_stats}", sep='\t')

       fraction_unmapped_reads = float(df_unmapped_reads_seqkit_stats['num_seqs'].iloc[0])/float(df_raw_reads_seqkit_stats['num_seqs'].iloc[0])
       n50_hifiasm_meta_assembly = df_hifiasm_meta_assembly_seqkit_stats['N50'].iloc[0]
       n50_metaflye_assembly = df_metaflye_assembly_seqkit_stats['N50'].iloc[0]
       n50_hicanu_assembly = df_hicanu_assembly_seqkit_stats['N50'].iloc[0]

       print(f"Stats for: {wildcards.env}")
       print(f"-----------------")
       print(f"fraction_unmapped_reads: {fraction_unmapped_reads}")
       print(f"n50_hifiasm_meta_assembly: {n50_hifiasm_meta_assembly}")
       print(f"n50_metaflye_assembly: {n50_metaflye_assembly}")
       print(f"n50_hicanu_assembly: {n50_hicanu_assembly}")