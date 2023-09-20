# latent-rna
 Extract latent transcriptomic phenotypes

## Usage

### 1. Define gene feature bins

Gene annotations should be preprocessed to merge the exons for all isoforms of a gene into one set of non-overlapping exon regions. To do this, use the [`collapse_annotation.py` script from the GTEx pipeline](https://github.com/broadinstitute/gtex-pipeline/blob/master/gene_model/collapse_annotation.py):

```shell
python scripts/collapse_annotation.py Homo_sapiens.GRCh38.106.gtf.gz Homo_sapiens.GRCh38.106.genes.gtf --collapse_only
gzip Homo_sapiens.GRCh38.106.genes.gtf
```

Then, split exons, introns, upstream, and downstream regions into bins:

```shell
python scripts/get_gene_bins.py --gtf Homo_sapiens.GRCh38.106.genes.gtf.gz --chromosomes chr_lengths.genome --output gene_bins.bed.gz
```

This step also filters genes to include only those with `gene_biotype`/`gene_type` of  "protein_coding", and excludes bins that overlap the exon of another gene. `chr_lengths.genome` is a table of chromosome names and lengths, e.g. chrNameLength.txt from the STAR index, to sort chromosomes.

#### Command line options

```
usage: get_gene_bins.py [-h] -g GTF -c CHROMOSOMES [-n N_BINS] -o OUTPUT

Get gene feature bins from GTF file

options:
  -h, --help            show this help message and exit
  -g GTF, --gtf GTF     GTF file
  -c CHROMOSOMES, --chromosomes CHROMOSOMES
                        Chromosome lengths file, e.g. chrNameLength.txt from STAR index, to sort chromosomes.
  -n N_BINS, --n-bins N_BINS
                        Number of bins to split each feature (exon, intron, etc.) into
  -o OUTPUT, --output OUTPUT
                        Output file (BED)
```

### 2. RNA-seq coverage counts

Run [`bedtools coverage` (`coverageBed`)](https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html) to count the reads overlapping each bin, e.g.:

```shell
bedtools coverage -split -sorted -counts -a ref/gene_bins.bed.gz -b /path/to/bams/sampleA.bam -g chr_lengths.genome | cut -f7 > covg/sampleA.txt
```

This produces a files for each sample containing the bin counts, one per line, corresponding to the lines of `gene_bins.bed.gz`. In subsequent steps, a list of these files will be loaded from a file, so generate that list:

```shell
awk '{print "covg/"$1".txt"}' samples.txt > covgfiles.txt
```

### 3. Generating latent phenotypes

There are three stages to generate latent phenotypes from coverage counts:

1. `prepare`: Per-sample coverage files are first assembled into a bin by sample matrix, split into batches of N genes each, and saved.
2. `fit`: Coverage count matrices are loaded, one batch at a time, and used to fit a PCA model for each gene. If a project involves multiple datasets, e.g. tissues, the coverage count matrices for each dataset can be loaded and concatenated, and the models fit on the concatenated data.
3. `transform`: Coverage count matrices are loaded, one batch at a time, and used to transform the data using the fitted PCA models. The transformed data is saved as a table of samples by phenotypes, i.e. multiple PCs per gene. If a project involves multiple datasets, the coverage count matrices for each dataset can be loaded and transformed separately using the same set of models, so that the phenotypes correspond across datasets.

`examples/` contains a bash script with examples of running these steps in different scenarios. It also contains a `Snakefile` that can be modified as needed and used to run the steps in a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow.

#### Command line options

```
usage: latent_RNA.py [-h] {prepare,fit,transform} ...

Fit and/or apply PCA model on feature bin coverage data

options:
  -h, --help            show this help message and exit

subcommands:
  {prepare,fit,transform}
                        Choose a subcommand
    prepare             Prepare batched input coverage files. Coverage counts will be assembled and
                        saved in batches, `fit` will be run on each batch separately, and
                        `transform` will run all batches and produce the combined output.
    fit                 Fit PCA model
    transform           Apply PCA transformation
```

```
usage: latent_RNA.py prepare [-h] -i FILE -r FILE -d DIR [--batch-size N]

options:
  -h, --help            show this help message and exit
  -i FILE, --inputs FILE
                        File containing paths to all per-sample coverage files. Each one has one
                        integer per line corresponding to rows of `regions`.
  -r FILE, --regions FILE
                        BED file containing regions to use for PCA. Must have start, end and region
                        ID in 2nd, 3rd, and 4th columns. Rows must correspond to rows of input
                        coverage files.
  -d DIR, --output-dir DIR
                        Directory where per-batch numpy binary files will be written.
  --batch-size N        Number of genes (at most) per batch. Default 200.
```

```
usage: latent_RNA.py fit [-h] (-d DIR [DIR ...] | --dir-file FILE) [-b N] -m DIR
                         [--n-samples-max N] [-v FLOAT] [-n N]

options:
  -h, --help            show this help message and exit
  -d DIR [DIR ...], --batch-covg-dir DIR [DIR ...]
                        Directory of per-batch numpy binary files. Specify multiple directories to
                        load data from all of them and fit models using the combined dataset. All
                        datasets must have been generated using the same gene bins file.
  --dir-file FILE       File containing list of directories of per-batch numpy binary files. Use
                        this instead of -d/--gene-covg-dir in case of many directories.
  -b N, --batch N       Batch number to load and fit. Batch numbers start from 0. If omitted, all
                        batches will be loaded and fit in sequence.
  -m DIR, --models-dir DIR
                        Directory to save PCA model pickle files
  --n-samples-max N     Max number of samples to use for fitting PCA models faster. Loaded samples
                        are randomly subsetted if higher than this, done after loading and
                        concatenating data from multiple datasets if applicable, and the sample
                        subsets are chosen independently per gene. Pass 0 for no cutoff.
  -v FLOAT, --var-expl-max FLOAT
                        Max variance explained by the PCs kept per gene. Pass 0 or 1 for no
                        variance explained cutoff. Default 0.8.
  -n N, --n-pcs-max N   Max number of PCs to keep per gene. Pass 0 for no cutoff. Default 32.
```

```
usage: latent_RNA.py transform [-h] -d DIR -m DIR -o FILE

options:
  -h, --help            show this help message and exit
  -d DIR, --batch-covg-dir DIR
                        Directory of per-batch numpy binary files.
  -m DIR, --models-dir DIR
                        Directory of saved PCA models (*.pickle) to load and use for transformation
  -o FILE, --output FILE
                        Output file (TSV)
```
