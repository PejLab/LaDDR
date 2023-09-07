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
