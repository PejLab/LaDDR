# latent-rna
 Extract latent transcriptomic phenotypes

## Usage

### Define gene feature bins

Gene annotations should be preprocessed to merge the exons for all isoforms of a gene into one set of non-overlapping exon regions. To do this, use the [`collapse_annotation.py` script from the GTEx pipeline](https://github.com/broadinstitute/gtex-pipeline/blob/master/gene_model/collapse_annotation.py):

```shell
python scripts/collapse_annotation.py Homo_sapiens.GRCh38.106.gtf.gz Homo_sapiens.GRCh38.106.genes.gtf --collapse_only
gzip Homo_sapiens.GRCh38.106.genes.gtf
```

Then, split exons, introns, upstream, and downstream regions into bins:

```shell
python scripts/get_gene_bins.py --gtf Homo_sapiens.GRCh38.106.genes.gtf.gz --chromosomes chr_lengths.genome --output gene_bins.bed.gz
```

This step also filters genes to include only those with `gene_biotype`/`gene_type` of  "protein_coding". `chr_lengths.genome` is a table of chromosome names and lengths, e.g. chrNameLength.txt from the STAR index, to sort chromosomes.

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

### RNA-seq coverage counts

Run [`bedtools coverage` (`coverageBed`)](https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html) to count the reads overlapping each bin, e.g.:

```shell
bedtools coverage -split -sorted -counts -a ref/gene_bins.bed.gz -b /path/to/bams/sampleA.bam -g chr_lengths.genome | cut -f7 > covg/sampleA.txt
```

This produces a files for each sample containing the bin counts, one per line, corresponding to the lines of `gene_bins.bed.gz`. In subsequent steps, a list of these files will be loaded from a file, so generate that list:

```shell
awk '{print "covg/"$1".txt"}' samples.txt > covgfiles.txt
```

### Generating latent phenotypes

Depending on the size of the dataset and usage needs, this same method can be run in different ways. First:

- The per-gene PCA models can be computed and saved, then applied to the same or other datasets. Applying the same set of models to multiple datasets allows the features to correspond across datasets.
- Each gene can be fit and transformed in one step, so that the input data is only loaded once and models aren't saved.

Second:

- For small datasets, all per-sample coverage files can be loaded directly, assembled and split by gene, and processed in sequence.
- For datasets where the coverage counts can't all be loaded into memory at once, per-sample coverage files are first assembled in batches into per-gene tables and saved. Coverage counts can then be loaded and processed in sequence per gene.
  <!-- - Or, the program can be run for a single gene, allowing genes to be processed in parallel using separate calls to the program. -->

Thus, the four ways to run the procedure are:

#### Small dataset, save and/or reuse models

```shell
python latent_RNA.py fit -i covgfiles.txt -r gene_bins.bed.gz -o models.pkl
python latent_RNA.py transform -i covgfiles.txt -r gene_bins.bed.gz -m models.pkl -o phenos.tsv.gz
```

#### Small dataset, one-time transformation

```shell
python latent_RNA.py fit-transform -i covgfiles.txt -r gene_bins.bed.gz -o phenos.tsv.gz
```

#### Large dataset, save and/or reuse models

```shell
python latent_RNA.py prepare -i covgfiles.txt -r gene_bins.bed.gz -d gene_covg/
python latent_RNA.py fit -d gene_covg/ -r gene_bins.bed.gz -o models.pkl
python latent_RNA.py transform -d gene_covg/ -r gene_bins.bed.gz -m models.pkl -o phenos.tsv.gz
```

#### Large dataset, one-time transformation

```shell
python latent_RNA.py prepare -i covgfiles.txt -r gene_bins.bed.gz -d gene_covg/
python latent_RNA.py fit-transform -d gene_covg/ -r gene_bins.bed.gz -o phenos.tsv.gz
```

#### Multiple datasets

To process multiple datasets, e.g. GTEx tissues, in a way that results in phenotypes that correspond across datasets, it is recommended to run the `prepare` step separately for each tissue, so that each dataset's per-gene coverage files can be used to run the `transform` step separately. The `fit` step can take a file listing multiple gene coverage directories, and will concatenate per gene to fit the models:

```shell
python latent_RNA.py prepare -i covgfiles.tissue1.txt -r gene_bins.bed.gz -d gene_covg_tissue1/
python latent_RNA.py prepare -i covgfiles.tissue2.txt -r gene_bins.bed.gz -d gene_covg_tissue2/
# Gene coverage directories can either be listed with -d, or if there are many, listed in a file:
echo 'gene_covg_tissue1/\ngene_covg_tissue2/' > gene_covg_dirs.txt
python latent_RNA.py fit -d gene_covg_tissue1/ gene_covg_tissue2/ -r gene_bins.bed.gz -o models.pkl
python latent_RNA.py transform -d gene_covg_tissue1/ -r gene_bins.bed.gz -m models.pkl -o phenos.tissue1.tsv.gz
python latent_RNA.py transform -d gene_covg_tissue2/ -r gene_bins.bed.gz -m models.pkl -o phenos.tissue2.tsv.gz
```

If using many tissues, you can use loops and pass a file listing the coverage directories to the fit step:

```shell
cat tissues.txt | while read tissue; do
    python latent_RNA.py prepare -i covgfiles.${tissue}.txt -r gene_bins.bed.gz -d gene_covg_${tissue}/
done
awk '{print "gene_covg_"$1}' tissues.txt > gene_covg_dirs.txt
python latent_RNA.py fit --dir-file gene_covg_dirs.txt -r gene_bins.bed.gz -o models.pkl
cat tissues.txt | while read tissue; do
    python latent_RNA.py transform -d gene_covg_${tissue}/ -r gene_bins.bed.gz -m models.pkl -o phenos.${tissue}.tsv.gz
done
```

#### Command line options

```
usage: latent_RNA.py [-h] {prepare,fit,transform,fit-transform} ...

Fit and/or apply PCA model on feature bin coverage data

options:
  -h, --help            show this help message and exit

subcommands:
  {prepare,fit,transform,fit-transform}
                        Choose a subcommand
    prepare             Prepare per-gene input coverage files
    fit                 Fit PCA model
    transform           Apply PCA transformation
    fit-transform       Fit PCA model and apply transformation without saving models
```

```
usage: latent_RNA.py prepare [-h] -i FILE -r FILE -d DIR

options:
  -h, --help            show this help message and exit
  -i FILE, --inputs FILE
                        File containing paths to all input coverage files.
  -r FILE, --regions FILE
                        BED file containing regions to use for PCA. Must have start, end and region ID in 2nd, 3rd, and 4th columns. Rows must correspond to rows of input coverage files.
  -d DIR, --output-dir DIR
                        Directory where per-gene numpy binary files will be written.
```

```
usage: latent_RNA.py fit [-h] (-i FILE | -d DIR [DIR ...] | --dir-file FILE) -r FILE [--n-samples-max N] [-v FLOAT] [-n N] -o FILE

options:
  -h, --help            show this help message and exit
  -i FILE, --inputs FILE
                        File containing paths to all input coverage files.
  -d DIR [DIR ...], --gene-covg-dir DIR [DIR ...]
                        Directory of per-gene numpy binary files. Specify multiple directories to load them all and treat as a single dataset.
  --dir-file FILE       File containing list of directories of per-gene numpy binary files. Use this instead of -d/--gene-covg-dir in case of many directories.
  -r FILE, --regions FILE
                        BED file containing regions to use for PCA. Must have start, end and region ID in 2nd, 3rd, and 4th columns. Rows must correspond to rows of input coverage files.
  --n-samples-max N     Max number of samples to use for fitting PCA models for efficiency. Used only for preprocessed per-gene inputs. Loaded samples are randomly subsetted if higher than this, done after
                        loading and concatenating data from multiple datasets if applicable, and the sample subsets are chosen independently per gene. Pass 0 for no cutoff. Default 1024.
  -v FLOAT, --var-expl-max FLOAT
                        Max variance explained by the PCs kept per gene. Pass 0 or 1 for no variance explained cutoff. Default 0.8.
  -n N, --n-pcs-max N   Max number of PCs to keep per gene. Pass 0 for no cutoff. Default 32.
  -o FILE, --output FILE
                        Output file (*.pickle) to save PCA models
```

```
usage: latent_RNA.py transform [-h] (-i FILE | -d DIR) -r FILE -m FILE -o FILE

options:
  -h, --help            show this help message and exit
  -i FILE, --inputs FILE
                        File containing paths to all input coverage files. Base file name before first "." is sample ID.
  -d DIR, --gene-covg-dir DIR
                        Directory of per-gene numpy binary files.
  -r FILE, --regions FILE
                        BED file containing regions to use for PCA. Must have start, end and region ID in 2nd, 3rd, and 4th columns. Rows must correspond to rows of input coverage files.
  -m FILE, --models FILE
                        PCA models file (*.pickle) to load and use for transformation
  -o FILE, --output FILE
                        Output file (TSV)
```

```
usage: latent_RNA.py fit-transform [-h] (-i FILE | -d DIR) -r FILE [-v FLOAT] [-n N] -o FILE

options:
  -h, --help            show this help message and exit
  -i FILE, --inputs FILE
                        File containing paths to all input coverage files. Base file name before first "." is sample ID.
  -d DIR, --gene-covg-dir DIR
                        Directory of per-gene numpy binary files.
  -r FILE, --regions FILE
                        BED file containing regions to use for PCA. Must have start, end and region ID in 2nd, 3rd, and 4th columns. Rows must correspond to rows of input coverage files.
  -v FLOAT, --var-expl-max FLOAT
                        Max variance explained by the PCs kept per gene. Pass 0 or 1 for no variance explained cutoff. Default 0.8.
  -n N, --n-pcs-max N   Max number of PCs to keep per gene. Pass 0 for no cutoff. Default 32.
  -o FILE, --output FILE
                        Output file (TSV)
```
