## Clone the repository
```
git clone https://github.com/ay-amityadav/HiFiBGC_analyses
```
## Download dataset
Download [HiFiBGC_0.1.13_Run.tar.gz](https://zenodo.org/records/10874958/files/HiFiBGC_0.1.13_Run.tar.gz?download=1), uncompress it, and put it under folder `HiFiBGC_analyses`.

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

