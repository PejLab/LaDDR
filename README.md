# latent-rna
 Extract latent transcriptomic phenotypes

## Usage

This method involves four main steps as outlined below. If using existing genomic bin definitions, only steps 2-4 are required. If also using pretrained latent phenotyping models, only steps 2 and 4 are required. `examples/` contains a bash script with examples of running these steps in different scenarios. It also contains a `Snakefile` that can be modified as needed and used to run the steps in a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow. We recommend using Snakemake and adapting this `Snakefile` for your project.

### 0. Get base-level RNA-seq coverage

This tool uses base-level RNA-seq coverage in bigWig format. For BAM files, first convert them to bigWig using [`bamCoverage`](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) from [`deeptools`](https://deeptools.readthedocs.io/en/develop/index.html):

```shell
bamCoverage -b /path/to/bams/sampleA.bam -o bigWig/sampleA.bw -of bigwig --binSize 1 -p 4
```

### 1. Define genomic bins

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
python scripts/get_gene_bins.py \
    --gtf Homo_sapiens.GRCh38.106.genes.gtf.gz \
    --chromosomes chr_lengths.genome \
    --outdir gene_bins/
```

This step also filters genes to include only those with `gene_biotype`/`gene_type` of  "protein_coding", and excludes bins that overlap the exon of another gene. `chr_lengths.genome` is a table of chromosome names and lengths, e.g. chrNameLength.txt from the STAR index, to sort chromosomes. Genes are assigned to batches so that the coverage processing and model fitting steps can be run in batches for parallelization and to reduce memory. In the final phenotyping step, all batches are recombined into one phenotype file.

#### Adaptive binning

The above binning method relies on exon definitions in the GTF file. Binning can also be determined using coverage data, partitioning each gene in a way that aims to define more, smaller bins in areas of greater variation across samples:

```shell
python scripts/get_gene_bins.py \
    --gtf Homo_sapiens.GRCh38.106.genes.gtf.gz \
    --chromosomes chr_lengths.genome \
    --binning-method adaptive3 \
    --bigwig-paths-file covg_bigwig_files.txt \
    --bins-per-gene 256 \
    --outdir gene_bins/
```

#### Command line options

```
usage: get_gene_bins.py [-h] -g FILE -c FILE --outdir PATH
                        [--binning-method {adaptive1,adaptive2,adaptive3,bin-width,n-bins}]
                        [--bigwig-paths-file FILE]
                        [--min-mean-total-covg FLOAT] [--max-corr FLOAT]
                        [--bins-per-gene N] [--bin-width-coding N]
                        [--bin-width-noncoding N] [--n-bins N]
                        [--max-bin-width N] [--batch-size N] [--batch N]

Partition genes into bins for summarizing coverage data

options:
  -h, --help            show this help message and exit
  -g FILE, --gtf FILE   Transcript annotation in GTF format. Must be the
                        collapsed annotation produced by
                        `collapse_annotation.py`.
  -c FILE, --chromosomes FILE
                        Chromosome lengths file, e.g. chrNameLength.txt from
                        STAR index, to sort chromosomes.
  --outdir PATH         Directory in which to save per-batch BED files.
  --binning-method {adaptive1,adaptive2,adaptive3,bin-width,n-bins}
                        Whether to determine bins adaptively from coverage
                        data, split all coding/noncoding regions into fixed-
                        width bins, or split each region into a fixed number
                        of bins. (default: n-bins)
  --bigwig-paths-file FILE
                        For "adaptive" methods, file containing list of paths
                        to per-sample bigWig files to use for coverage data.
  --min-mean-total-covg FLOAT
                        For method "adaptive1", minimum allowed mean total
                        coverage per sample for a bin. (default: 128)
  --max-corr FLOAT      For method "adaptive1", maximum allowed correlation
                        between normalized coverage of adjacent bins.
                        (default: 0.8)
  --bins-per-gene N     For method "adaptive2" or "adaptive3", approximate
                        number of bins to create per gene on average.
                        (default: 128)
  --bin-width-coding N  For method "bin-width", width of bins for coding
                        (exonic) regions of the genes. (default: 16)
  --bin-width-noncoding N
                        For method "bin-width", width of bins for noncoding
                        (non-exonic) regions of the genes. (default: 128)
  --n-bins N            For method "n-bins", number of bins to split each
                        feature (exon, intron, etc.) into. (default: 24)
  --max-bin-width N     After a binning method is run, any bins larger than
                        this will be split up. (default: 1024)
  --batch-size N        Number of genes (at most) per batch. (default: 200)
  --batch N             Batch ID to process. Batch IDs are integers starting
                        from 0. If omitted, all batches will be processed.
```

### 2. Bin and normalize RNA-seq coverage


For each gene batch, use its bin definitions to get mean coverage per bin for all samples and normalize:

```shell
python latent_RNA.py prepare -i covg_bigwig_files.txt --bins-dir gene_bins/ -b 0 -o covg_norm/
```

At this stage you can also provide quantified explicit phenotypes, e.g. from [Pantry](https://github.com/PejLab/Pantry), to regress out. Training and applying models on this residualized coverage data results in "residual" latent RNA phenotypes, which can complement the explicit phenotypes by representing uncharacterized transcriptomic variation.

#### Command line options

```
usage: latent_RNA.py prepare [-h] -i FILE --bins-dir DIR
                             (-b N | --n-batches N) [-p FILE [FILE ...] |
                             --pheno-paths-file FILE] -o DIR

options:
  -h, --help            show this help message and exit
  -i FILE, --bigwig-paths-file FILE
                        File containing list of paths to per-sample bigWig
                        files. Basenames of files will be used as sample IDs.
  --bins-dir DIR        Directory of per-batch BED files containing bin
                        regions. Must have start, end and bin ID in 2nd, 3rd,
                        and 4th columns. Rows must correspond to rows of
                        corresponding input coverage file.
  -b N, --batch N       Batch ID to process. Batch IDs are integers starting
                        from 0.
  --n-batches N         To load and fit all batches in sequence, provide
                        number of batches instead of a specific batch.
  -p FILE [FILE ...], --pheno-paths FILE [FILE ...]
                        One or more paths to phenotype tables to regress out
                        of the coverage data per gene prior to model input.
                        Files should be in bed format, i.e. input format for
                        tensorqtl. Gene IDs are parsed from the 4th column
                        from the start up to the first non-alphanumeric
                        character.
  --pheno-paths-file FILE
                        File containing list of paths to phenotype tables.
  -o DIR, --output-dir DIR
                        Directory where per-batch numpy binary files with
                        normalized coverage will be written.
```

### 3. Fit latent RNA phenotype models

Normalized coverage data are loaded, one batch at a time, and used to fit an FPCA or PCA model for each gene. If a project involves multiple datasets, e.g. tissues, the coverage count matrices for each dataset can be loaded together, and the models will be fit on the concatenated data.

```shell
python latent_RNA.py fit -d covg_norm/ -b 0 -m models/
```

#### Command line options

```
usage: latent_RNA.py fit [-h] (-d DIR [DIR ...] | --norm-covg-dir-file FILE)
                         (-b N | --n-batches N) -m DIR [-v FLOAT] [-n N]
                         [--use-fpca] [--fpca-x-values STRING]
                         [--fpca-basis STRING]

options:
  -h, --help            show this help message and exit
  -d DIR [DIR ...], --norm-covg-dir DIR [DIR ...]
                        Directory of per-batch numpy binary files with
                        normalized coverage. Specify multiple directories to
                        load data from all of them and fit models using the
                        combined dataset. All datasets must have been
                        generated using the same per-batch gene bins files.
  --norm-covg-dir-file FILE
                        File containing list of directories of per-batch numpy
                        binary files. Use this instead of -d/--norm-covg-dir
                        in case of many directories.
  -b N, --batch N       Batch ID to load and fit. Batch IDs are integers
                        starting from 0.
  --n-batches N         To load and fit all batches in sequence, provide
                        number of batches instead of a specific batch.
  -m DIR, --models-dir DIR
                        Directory in which to save model pickle files.
  -v FLOAT, --var-expl-max FLOAT
                        Max variance explained by the PCs kept per gene. Pass
                        0 or 1 for no variance explained cutoff. (default:
                        0.8)
  -n N, --n-pcs-max N   Max number of PCs to keep per gene. Pass 0 for no
                        cutoff. (default: 16)
  --use-fpca            Use functional PCA instead of regular PCA.
  --fpca-x-values STRING
                        Whether to use bin numbers or genomic positions as
                        x-values for functional PCA. (default: bin)
  --fpca-basis STRING   Basis function to use for functional PCA. `discrete`
                        will run discretized FPCA directly on the data,
                        `spline` will use a 4th-order B-spline basis.
                        (default: discrete)
```

### 4. Generate latent RNA phenotypes

Normalized coverage matrices are loaded, one batch at a time, and transformed using the fitted models. The transformed data is saved as a table of phenotypes by samples, i.e. multiple PCs per gene. If a project involves multiple datasets, the coverage count matrices for each dataset can be loaded and transformed separately using the same set of models, so that the phenotypes correspond across datasets.

```shell
python latent_RNA.py transform -d covg_norm/ -m models/ --n-batches 3 -o latent_phenos.tsv.gz
```

#### Command line options

```
usage: latent_RNA.py transform [-h] -d DIR -m DIR --n-batches N -o FILE

options:
  -h, --help            show this help message and exit
  -d DIR, --norm-covg-dir DIR
                        Directory of per-batch numpy binary files with
                        normalized coverage.
  -m DIR, --models-dir DIR
                        Directory of saved models (*.pickle) to load and use
                        for transformation.
  --n-batches N         Number of batches in the data. Latent phenotypes from
                        all batches will be computed and concatenated.
  -o FILE, --output FILE
                        Output file (TSV).
```
