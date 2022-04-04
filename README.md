# Snakemake workflow: Predicting gene expression using epigenetic signals captured by fragmentation patterns of cell-free DNA and related features

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)

This repository contains snakemake workflows for predicting gene expression using epigenetic signals captured by fragmentation patterns of cell-free DNA and related features. 


## Authors

- Kristin Köhler (@koehlek99)

## Usage


### Step 1: Obtain a copy of this workflow

[Clone](https://help.github.com/en/articles/cloning-a-repository) the repository to your local system, into the place where you want to perform the data analysis.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config_downsampling.yaml` to configure the coverage the samples are subsampled to and `samples_downsampling.tsv` to specify your sample setup.

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

```bash
conda create -c bioconda -c conda-forge -n snakemake snakemake
```

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

The workflows are executed from the repository root folder. The different analyses have to be executed separately. To specify the respective workflow use the `-s` switch followed by the path of the `Snakefile` (`workflow/Snakefile.smk`)

Activate the conda environment:

```bash
conda activate snakemake
```

Test your configuration by performing a dry-run via

```bash
snakemake -s workflow/Snakefile.smk --use-conda -n
```

Execute the workflow locally via

```bash
snakemake -s workflow/Snakefile.smk --use-conda --cores $N
```


## Workflows

### Downsampling

The downsampling step is contained in the `downsampling.smk` workflow.

#### Description
After calculating the coverage of each BAM file using mosdepth, the BAM files are subsampled to the (in config_downsampling.yaml) specified coverage. The resulting BAM files are indexed afterwards using samtools.   

#### Input

- configured by the user ([samples_downsampling.tsv](config/samples_downsampling.tsv)):
    - analysis ID
    - samples
    - path to sample .bam files
    - reference samples for plotting
    - genome build per sample
- configured by the user ([config_downsampling.yaml](config/config_downsampling.yaml.tsv)):
    - path to sample.tsv 
    - coverage to subsample to (cannot be larger than minimal coverage of all samples)

#### Output

- coverage statistics generated with mosdepth (results/downsampling/)
- BAM files with specified coverage and corresponding index files (results/downsampling/{GENOME}/)

### Gene Expression Analysis

The analysis is contained in the `snakefile_GE_analysis.smk` workflow.

#### Description

This workflow was taken from the gene expression analysis workflow from https://github.com/kircherlab/cfDNA.

#### Input

- included in the repository:
    - annotations
    - labels
    - RNAtable from Protein Atlas
        - protein atlas tissues + cell lines ["Extended"]

#### Output

- fft_summary tables (results/intermediate/healthy/FFT/)
- (normalized) WPS tables (results/intermediate/healthy/table/)
- (normalized) COV tables (results/intermediate/healthy/table/)

### Random forest prediction

The analysis is contained in the `rf_prediction.smk` workflow.

#### Description

This workflow extracts features from the previous generated tables (FFT, WPS, COV) and other resources (GC content, mean expression across tissues). The resulting feature tables are introduced to the provided Random forest model (resources/rf/rf_model.joblib) that was trained on two healthy samples using the expression of monocytes as reference labels. 

#### Input

- FFT tables (results/intermediate/{ID}/FFT_table/transcriptanno-{SAMPLE}-FFT_table.{COV}x.{GENOME}.tsv.gz)
- WPS tables (results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_WPS_normalized.{COV}x.{GENOME}.tsv.gz)
- COV tables (results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_COV_normalized.{COV}x.{GENOME}.tsv.gz)
- GC content (resources/rf/GC_content.csv)
- mean expression across tissues (resources/rf/mean_expression.csv)
- monocyte expression (resources/rf/RNAtableExtended.tsv.gz)


#### Output

- predicted expression per gene (results/predictions/{ID}/{GENOME}/)