## Clone the repository
```
git clone https://github.com/ay-amityadav/HiFiBGC_analyses
```
## Download dataset
Download `HiFiBGC_0.1.13_Run.tar.gz` from [figshare](https://doi.org/10.6084/m9.figshare.27194043.v3) or [zenodo](https://zenodo.org/records/10874958), uncompress it, and put it under folder `HiFiBGC_analyses`.

## Setup conda environments
```
conda env create --file envs/snakemake.yml

conda env create --file envs/jupyterlab.yml
```

## Run

### Comparison between four methods
Run 
```
conda activate snakemake

snakemake -s comparison_between_four_methods.smk --use-conda --cores 8 -p
```
Run
```
conda activate jupyterlab

jupyter lab
```
Thereafter, execute `comparison_between_four_methods.ipynb`

### Clinker analysis
Run
```
conda activate snakemake

snakemake -s clinker.smk --use-conda --cores 8 -p
```
Run
```
conda activate jupyterlab

jupyter lab
```
Thereafter, execute `clinker_analysis.ipynb` with jupyter kernel created from conda environment `envs/raincloudplots.yml`.

### Comparison within HiFiBGC (Upsetplots) 
Run
```
conda activate snakemake

snakemake -s upsetplot.smk --use-conda --cores 8 -p
``` 

### Complete BGC from unmapped reads
Run
```
conda activate jupyterlab

jupyter lab
```
Thereafter, execute `complete_BGCs_from_unmapped-reads.ipynb` with jupyter kernel created from conda environment `envs/dna_features_viewer.yml`.

### Sequence statistics
Download raw-read files for SRA-Ids: SRR15275213 (Human), ERR7015089 (Sludge), SRR10963010 (Sheep) and SRR15214153 (Chicken). And put them under the folder `raw` with names `human.fastq`, `sludge.fastq`, `sheep.fastq` and `chicken.fastq`. 

Run
```
conda activate snakemake

snakemake -s sequence_stats.smk --use-conda --cores 8 -p
```

### BGC database statistics
Download [mibig_gbk_3.1.tar.gz](https://dl.secondarymetabolites.org/mibig/mibig_gbk_3.1.tar.gz), and uncompress it, and put it under folder `HiFiBGC_analyses`.

Download [BiG-SLICE database](https://s3.ap-northeast-1.wasabisys.com/gigadb-datasets/live/pub/10.5524/100001_101000/100826/data/full_run_result.zip), and uncompress it, and put it under folder `HiFiBGC_analyses`.

Run
```
conda activate jupyterlab

jupyter lab
```
Thereafter, execute `BGC_databases.ipynb`


