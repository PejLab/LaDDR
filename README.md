# latent-rna
 Extract latent transcriptomic phenotypes

## Usage

### RNA-seq coverage counts



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
python latent_RNA.py fit --dir-file gene_covg_dirs.txt -r gene_bins.bed.gz -o models.pkl
python latent_RNA.py transform -d gene_covg_tissue1/ -r gene_bins.bed.gz -m models.pkl -o phenos_tissue1.tsv.gz
python latent_RNA.py transform -d gene_covg_tissue2/ -r gene_bins.bed.gz -m models.pkl -o phenos_tissue2.tsv.gz
```
