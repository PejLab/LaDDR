# latent-rna
 Extract latent transcriptomic phenotypes

## Installation

```shell
cd latent-rna
pip install -e .
```

## Setup

Run this command to generate a project directory `myproject` with the basic config and to template options for running the pipeline:

```shell
latent-rna init myproject
```

```
Project initialized at myproject

Files created:
  myproject/
    config.yaml: Edit this file to configure your project
    run.sh: Shell script for running the project
    Snakefile: Snakemake file for running the project
```

Either `run.sh` or `Snakefile` can be modified as needed to run the project. We recommend adapting the `Snakefile` for your project and using the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow system.

The default config contains parameters for the recommended binning and model fitting methods. To start with a config file containing all method options and parameters, use `latent-rna init myproject --config-type extended`.

`examples/` contains two project directories set up to run example projects on two small datasets using either method. If you install the package from pip and do not have `examples/`, you can generate its two subdirectories with these commands:

```shell
latent-rna init shell --config-type example --template shell
latent-rna init snakemake --config-type example --template snakemake
```

### Get base-level RNA-seq coverage

This tool uses base-level RNA-seq coverage in bigWig format. For BAM files, first convert them to bigWig using [`bamCoverage`](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) from [`deeptools`](https://deeptools.readthedocs.io/en/develop/index.html):

```shell
bamCoverage -b /path/to/bams/sampleA.bam -o covg_bigwig/sampleA.bw -of bigwig --binSize 1 -p 4
```

Then, make a coverage manifest file listing the dataset, sample ID, and path to all bigWig files, e.g.:

```
tissue1	S000	../data_input/covg_bigwig/tissue1/S000.bw
tissue1	S001	../data_input/covg_bigwig/tissue1/S001.bw
...
tissue8	S098	../data_input/covg_bigwig/tissue8/S098.bw
tissue8	S099	../data_input/covg_bigwig/tissue8/S099.bw
```

A latent phenotype table will be produced for each dataset containing phenotypes for its samples. The reason multiple datasets are supported is that the coverage data can be loaded together and the latent phenotype models can be fit jointly, so that the models are trained on biologically diverse data and the phenotypes correspond across datasets.

## Usage

This method involves four main steps as outlined below. If using existing genomic bin definitions, only steps 2-4 are required. If also using pretrained latent phenotyping models, only steps 2 and 4 are required.

### 1. Define genomic bins

The binning step requires two reference files:

- Chromosome lengths, which can be `chrNameLength.txt` from the STAR index or the `*.fa.fai` genome FASTA index file.
- Gene annotations, preprocessed to merge the exons for all isoforms of a gene into one set of non-overlapping exon regions. To do this, use the [`collapse_annotation.py` script from the GTEx pipeline](https://github.com/broadinstitute/gtex-pipeline/blob/master/gene_model/collapse_annotation.py):

```shell
python scripts/collapse_annotation.py Homo_sapiens.GRCh38.106.gtf.gz collapsed.gtf --collapse_only
gzip collapsed.gtf
```

Then, split gene regions into bins to be used for summarizing and representing coverage data:

```shell
latent-rna binning --batch 0
```

This step also filters genes to include only those with `gene_biotype`/`gene_type` of  "protein_coding", and excludes bins that overlap the exon of another gene. Genes are assigned to batches so that the binning, coverage processing, and model fitting steps can be run in batches for parallelization and to reduce memory. In the final phenotyping step, all batches are recombined into one phenotype file. By default, binning is determined using coverage data, partitioning each gene in a way that aims to define more, smaller bins in areas of greater variation across samples.

#### Command line options

```
usage: latent-rna binning [-h] [-c FILE] [-p PROJECT_DIR] [-b N]

options:
  -h, --help            show this help message and exit
  -c FILE, --config FILE
                        Path to project configuration file. Defaults to
                        config.yaml in project directory.
  -p PROJECT_DIR, --project-dir PROJECT_DIR
                        Project directory. Paths in config are relative to
                        this. Defaults to current directory.
  -b N, --batch N       Batch ID to process. Batch IDs are integers starting
                        from 0. If omitted, all batches will be processed.
```

### 2. Bin and normalize RNA-seq coverage


For each gene batch, use its bin definitions to get mean coverage per bin for all samples and normalize:

```shell
latent-rna prepare --dataset dset1 --batch 0
```

At this stage you can also provide quantified explicit phenotypes, e.g. from [Pantry](https://github.com/PejLab/Pantry), to regress out. Training and applying models on this residualized coverage data results in "residual" latent RNA phenotypes, which can complement the explicit phenotypes by representing uncharacterized transcriptomic variation.

#### Command line options

```
usage: latent-rna prepare [-h] [-c FILE] [-p PROJECT_DIR] [-d NAME] [-b N]

options:
  -h, --help            show this help message and exit
  -c FILE, --config FILE
                        Path to project configuration file. Defaults to
                        config.yaml in project directory.
  -p PROJECT_DIR, --project-dir PROJECT_DIR
                        Project directory. Paths in config are relative to
                        this. Defaults to current directory.
  -d NAME, --dataset NAME
                        Name of dataset to process.
  -b N, --batch N       Batch ID to process. Batch IDs are integers starting
                        from 0. If omitted, all batches will be processed.
```

### 3. Fit latent RNA phenotype models

Normalized coverage data are loaded, one batch at a time, and used to fit an FPCA or PCA model for each gene. If a project involves multiple datasets, e.g. tissues, the coverage count matrices for each dataset can be loaded together, and the models will be fit on the concatenated data.

```shell
latent-rna fit --batch 0
```

#### Command line options

```
usage: latent-rna fit [-h] [-c FILE] [-p PROJECT_DIR] [-b N]

options:
  -h, --help            show this help message and exit
  -c FILE, --config FILE
                        Path to project configuration file. Defaults to
                        config.yaml in project directory.
  -p PROJECT_DIR, --project-dir PROJECT_DIR
                        Project directory. Paths in config are relative to
                        this. Defaults to current directory.
  -b N, --batch N       Batch ID to process. Batch IDs are integers starting
                        from 0. If omitted, all batches will be processed.
```

### 4. Generate latent RNA phenotypes

Normalized coverage matrices are loaded, one batch at a time, and transformed using the fitted models. The transformed data is saved as a table of phenotypes by samples, i.e. multiple PCs per gene. If a project involves multiple datasets, the coverage count matrices for each dataset can be loaded and transformed separately using the same set of models, so that the phenotypes correspond across datasets.

```shell
latent-rna transform --dataset dset1
```

#### Command line options

```
usage: latent-rna transform [-h] [-c FILE] [-p PROJECT_DIR] [-d NAME]

options:
  -h, --help            show this help message and exit
  -c FILE, --config FILE
                        Path to project configuration file. Defaults to
                        config.yaml in project directory.
  -p PROJECT_DIR, --project-dir PROJECT_DIR
                        Project directory. Paths in config are relative to
                        this. Defaults to current directory.
  -d NAME, --dataset NAME
                        Name of dataset to process.
```
