# latent-rna
 Extract latent transcriptomic phenotypes

## Installation

```shell
cd latent-rna
pip install -e .
```

## Usage

This method involves four main steps as outlined below. If using existing genomic bin definitions, only steps 2-4 are required. If also using pretrained latent phenotyping models, only steps 2 and 4 are required. `examples/` contains a bash script with examples of running these steps in different scenarios. It also contains a `Snakefile` that can be modified as needed and used to run the steps in a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow. We recommend using Snakemake and adapting this `Snakefile` for your project.

### 0. Get base-level RNA-seq coverage

This tool uses base-level RNA-seq coverage in bigWig format. For BAM files, first convert them to bigWig using [`bamCoverage`](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) from [`deeptools`](https://deeptools.readthedocs.io/en/develop/index.html):

```shell
bamCoverage -b /path/to/bams/sampleA.bam -o bigWig/sampleA.bw -of bigwig --binSize 1 -p 4
```

### 1. Define genomic bins

Chromosome lengths are needed for the binning step, and can be extracted from the genome FASTA index file:

```shell
cut -f1,2 Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai > chr_lengths.genome
```

Gene annotations should be preprocessed to merge the exons for all isoforms of a gene into one set of non-overlapping exon regions. To do this, use the [`collapse_annotation.py` script from the GTEx pipeline](https://github.com/broadinstitute/gtex-pipeline/blob/master/gene_model/collapse_annotation.py):

```shell
python scripts/collapse_annotation.py \
    Homo_sapiens.GRCh38.106.gtf.gz \
    Homo_sapiens.GRCh38.106.genes.gtf \
    --collapse_only
gzip Homo_sapiens.GRCh38.106.genes.gtf
```

Then, split exons, introns, upstream, and downstream regions into bins:

```shell
latent-rna binning \
    --gtf Homo_sapiens.GRCh38.106.genes.gtf.gz \
    --chromosomes chr_lengths.genome \
    --outdir gene_bins/
```

This step also filters genes to include only those with `gene_biotype`/`gene_type` of  "protein_coding", and excludes bins that overlap the exon of another gene. `chr_lengths.genome` is a table of chromosome names and lengths, e.g. chrNameLength.txt from the STAR index, to sort chromosomes. Genes are assigned to batches so that the coverage processing and model fitting steps can be run in batches for parallelization and to reduce memory. In the final phenotyping step, all batches are recombined into one phenotype file.

#### Adaptive binning

The above binning method relies on exon definitions in the GTF file. Binning can also be determined using coverage data, partitioning each gene in a way that aims to define more, smaller bins in areas of greater variation across samples:

```shell
latent-rna binning --batch 0
```

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
