import numpy as np
import pandas as pd
from typing import Iterator
from pathlib import Path
import pyBigWig
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm

class CoverageData:
    def __init__(self, norm_covg_dirs: list, batch_id: int):
        """Load normalized coverage data for one batch

        Args:
            norm_covg_dirs: List of directories containing normalized coverage
              data, e.g. to train models on multiple datasets. To load a single
              dataset, pass a list with one directory.
            batch_id: ID of the batch to load.
        """
        self.samples = []
        for d in norm_covg_dirs:
            with open(d / 'samples.txt') as f:
                self.samples.extend([l.strip() for l in f])
        if len(norm_covg_dirs) > 1:
            # Make sample IDs unique for concatenation, and ID won't be saved in the models anyway.
            self.samples = [f'S{i + 1}' for i in range(len(self.samples))]
        self.coverage, self.bins = self.load_coverage(norm_covg_dirs, batch_id)
        self.genes = list(self.bins.groupby('gene_id').groups.keys())

    def load_coverage(self, norm_covg_dirs: list, batch_id: int) -> tuple:
        """Load normalized coverage data for one batch
        
        Args:
            norm_covg_dirs: List of directories containing normalized coverage
              data.
            batch_id: ID of the batch to load.

        Returns:
            Tuple of coverage DataFrame and bins DataFrame.
        """
        mats = [np.load(d / f'batch_{batch_id}.npy') for d in norm_covg_dirs]
        mat = np.concatenate(mats, axis=1)
        bin_file = norm_covg_dirs[0] / f'batch_{batch_id}.bins.tsv.gz'
        bins = pd.read_csv(bin_file, sep='\t', index_col=[0, 1])
        assert mat.shape[0] == bins.shape[0]
        assert mat.shape[1] == len(self.samples)
        df = pd.DataFrame(mat, index=bins.index, columns=self.samples)
        return df, bins

    def by_gene(self) -> Iterator[tuple]:
        """Iterate over coverage data for each gene
        
        Yields:
            Tuple of gene ID and coverage DataFrame for the gene
        """
        for gene_id, df in self.coverage.groupby('gene_id'):
            yield gene_id, df

def load_bins(binfile: Path) -> pd.DataFrame:
    """Load bin information from a BED file

    The genomic coordinates will be used to extract mean coverage per bin from
    bigWig files. The gene-relative positions encoded in the bin IDs will be
    used for normalization and modeling.

    Args:
        binfile: Path to the BED file containing bin information. The 4th column
          (bin ID) will be parsed, and must contain gene ID, start (relative to
          gene start), and end separated by underscores.

    Returns:
        The DataFrame containing bin information. Rows are bins and are indexed
        by gene_id and pos (of bin center, relative to gene start), and columns
        are chrom, chrom_start, and chrom_end.
    """
    bins = pd.read_csv(
        binfile,
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 3],
        names=["chrom", "chrom_start", "chrom_end", "bin"],
        index_col="bin",
    )
    bins.index = bins.index.str.split("_", n=2, expand=True)
    bins.index.set_names(['gene_id', 'start', 'end'], inplace=True)
    bins = bins.reset_index()
    bins['chrom'] = bins['chrom'].astype(str)
    bins['start'] = bins['start'].astype(int)
    bins['end'] = bins['end'].astype(int)
    bins['length'] = bins['end'] - bins['start']
    bins['pos'] = bins['start'] + bins['length'] // 2
    bins = bins.set_index(['gene_id', 'pos'])
    bins = bins[['chrom', 'chrom_start', 'chrom_end']]
    return bins

def covg_from_bigwigs(bigwig_manifest: pd.DataFrame, bins: pd.DataFrame) -> pd.DataFrame:
    """Load coverage data from bigWig files
    
    Args:
        bigwig_manifest: DataFrame containing bigWig manifest. Must have columns
          sample and path.
        bins: DataFrame containing bin information. Rows are bins, indexed by
          gene_id and pos, and columns include chrom, chrom_start, and
          chrom_end.

    Returns:
        The coverage DataFrame. Rows are bins and columns are sample IDs.
    """
    covg = np.zeros((bins.shape[0], len(bigwig_manifest)))
    for i, (_, row) in enumerate(bigwig_manifest.iterrows()):
        with pyBigWig.open(str(row['path'])) as bw:
            for j, (gene_id, pos) in enumerate(bins.index):
                chrom = str(bins.loc[(gene_id, pos), 'chrom'])
                start = bins.loc[(gene_id, pos), 'chrom_start']
                end = bins.loc[(gene_id, pos), 'chrom_end']
                try:
                    covg[j, i] = bw.stats(chrom, start, end, type='mean', exact=True)[0]
                except RuntimeError as e:
                    print(f"RuntimeError for {gene_id}, {chrom}:{start}-{end} in file {row['path']}: {e}", flush=True)
    samples = bigwig_manifest['sample'].tolist()
    df = pd.DataFrame(covg, index=bins.index, columns=samples)
    return df

def normalize_coverage(
        df: pd.DataFrame,
        bins: pd.DataFrame,
        scaling_factors: pd.Series,
        pseudocount: float = 8
) -> pd.DataFrame:
    """Normalize coverage for each gene in one batch.

    1. Adjust for sequencing depth:
       - Mean coverage per base per bin is converted to total coverage per bin.
       - Coverage is scaled by pre-computed sample-specific factors.
       - Coverage is normalized back to mean coverage per base for each bin.
    2. A pseudocount is added to allow log-transformation and reduce the impact
       of noise in low-coverage bins.
    3. Coverage values are log-transformed for variance stabilization.

    Args:
        df: The input DataFrame containing mean coverage per bin per sample.
          Rows are bins and columns are samples. Multi-indexed by gene_id and
          pos.
        bins: The input DataFrame containing bin information.
        scaling_factors: A Series containing scaling factors indexed by sample
          IDs.
        pseudocount: Pseudocount to add to coverage values (mean coverage per
          bp for each bin) before normalization.

    Returns:
        The normalized coverage DataFrame.

    Raises:
        AssertionError: If the input DataFrame contains any NaN values.
    """
    assert df.index.equals(bins.index)
    # Convert to total coverage per bin
    lengths = bins['chrom_end'] - bins['chrom_start']
    df = df.mul(lengths, axis=0)
    # Scale by pre-computed sample-specific factors
    df = df.div(scaling_factors[df.columns], axis=1)
    assert not df.isna().any().any()
    # Normalize back to mean coverage per bp for each bin
    df = df.div(lengths, axis=0)
    # Enable log-transformation and reduce impact of noise
    df += pseudocount
    # Variance stabilization
    df = np.log2(df)
    return df

def load_phenotypes(pheno_file: Path, samples: list) -> pd.DataFrame:
    """Load phenotype table in BED format
    
    Gene IDs are parsed from the 4th column from the start up to the first
    non-alphanumeric character. Chromosome, start, and end columns are ignored.

    Args:
        pheno_file: Path to the phenotype table file in BED format.
        samples: List of sample IDs.

    Returns:
        The phenotype table DataFrame. Rows are phenotypes and columns are
        samples in the same order as the input list.
    
    Raises:
        AssertionError: If the phenotype table is missing any samples
    """
    df = pd.read_csv(pheno_file, sep='\t', index_col=3, dtype={0: str, 1: str, 2: str})
    df = df.drop(df.columns[:3], axis=1)
    df.index.name = 'phenotype_id'
    # Check that all samples are present:
    missing = set(samples) - set(df.columns)
    assert len(missing) == 0, f'Phenotype table is missing samples: {missing}'
    df = df[samples]
    return df

def regress_out_phenos_single(y: pd.Series, x: np.array) -> pd.DataFrame:
    """Regress out phenotypes from one feature (bin)

    Run regression in two stages:
    1. Identify nominally significant regression features (p < 0.01) using all
      phenotypes.
    2. Run regression using only the nominally significant features.
    
    Args:
        y: The input Series containing coverage values for one bin.
        x: The input array containing phenotypes for the gene.

    Returns:
        The residuals of the regression.
    """
    # Subset to nominally significant regression features
    x_const = sm.add_constant(x)
    model_all = sm.OLS(y, x_const).fit()
    assert model_all.pvalues.index[0] == 'const'
    x_sig = x[:, model_all.pvalues[1:] < 0.01] # exclude intercept
    if x_sig.shape[1] == 0:
        return y
    model = LinearRegression().fit(x_sig, y)
    resid = y - model.predict(x_sig)
    return pd.Series(resid, index=y.index)

def regress_out_phenos_gene(df: pd.DataFrame, phenos: pd.core.groupby.DataFrameGroupBy, gene_id: str) -> pd.DataFrame:
    """Regress out phenotypes per gene
    
    Args:
        df: The input DataFrame containing normalized coverage per bin per
          sample for one gene.
        phenos: The input DataFrameGroupBy containing phenotypes to regress out.
        gene_id: The gene ID to process.

    Returns:
        The DataFrame of residuals after regressing out phenotypes.
    """
    if gene_id in phenos.groups:
        x = phenos.get_group(gene_id).to_numpy().T
        df = df.apply(regress_out_phenos_single, axis=1, args=(x,)) # axis=1 means by row, i.e. index of Series is column names
    return df

def regress_out_phenos(covg: pd.DataFrame, pheno_files: list) -> pd.DataFrame:
    """Regress out phenotypes from coverage data
    
    Args:
        covg: The input DataFrame containing normalized coverage per bin per
          sample for one gene.
        pheno_files: List of paths to phenotype tables.

    Returns:
        The DataFrame of residuals after regressing out phenotypes.
    """
    print(f'Loading phenotypes from {len(pheno_files)} file{"" if len(pheno_files) == 1 else "s"} to regress out...', flush=True)
    phenos = [load_phenotypes(f, covg.columns).reset_index() for f in pheno_files]
    phenos = pd.concat(phenos, axis=0)
    # Parse gene IDs from phenotype IDs:
    pheno_genes = phenos['phenotype_id'].str.extract(r'([a-zA-Z0-9]+)')[0]
    phenos = phenos.drop('phenotype_id', axis=1)
    phenos = phenos.groupby(pheno_genes)
    # Regress out phenotypes:
    prev_index = covg.index.copy()
    covg = covg.groupby("gene_id", group_keys=False).apply(
        lambda c: regress_out_phenos_gene(c, phenos, c.name)
    )
    assert covg.index.equals(prev_index)
    return covg

def prepare_coverage(
        bigwig_manifest: pd.DataFrame,
        bins_dir: Path,
        scaling_factors: pd.Series,
        batches: list,
        pheno_files: list,
        outdir: Path
) -> None:
    """Assemble per-sample coverage, split into batches, and save to numpy binary files
    
    Args:
        bigwig_manifest: DataFrame containing bigWig manifest. Must have columns
          sample and path.
        bins_dir: Directory of per-batch BED files containing bin regions. Must
          have start, end and bin ID in 2nd, 3rd, and 4th columns. Rows must
          correspond to rows of corresponding input coverage file.
        scaling_factors: A Series containing scaling factors indexed by sample
          IDs.
        batches: List of batch IDs to process. Batch IDs are integers starting
          from 0.
        pheno_files: List of paths to phenotype tables to regress out of the
          coverage data per gene prior to model input.
        outdir: Directory where per-batch numpy binary files with normalized
          coverage will be written.
    """
    for batch in batches:
        print(f'=== Preparing coverage for batch {batch} ===', flush=True)
        binfile = bins_dir / f'batch_{batch}.bed.gz'
        bins = load_bins(binfile)
        covg = covg_from_bigwigs(bigwig_manifest, bins)
        # Sort by gene and position along gene instead of genomic coordinate
        covg = covg.sort_index()
        bins = bins.loc[covg.index, :]
        covg = normalize_coverage(covg, bins, scaling_factors)
        if len(pheno_files) > 0:
            covg = regress_out_phenos(covg, pheno_files)
        outdir.mkdir(exist_ok=True)
        mat = covg.to_numpy()
        np.save(outdir / f'batch_{batch}.npy', mat)
        b = bins.loc[covg.index, :]
        b.to_csv(outdir / f'batch_{batch}.bins.tsv.gz', sep='\t')
        # If samples.txt doesn't exist, write it:
        if not (outdir / 'samples.txt').exists():
            with open(outdir / 'samples.txt', 'w') as f:
                f.write('\n'.join(covg.columns) + '\n')
        print(f'Coverage saved in {outdir}', flush=True)

