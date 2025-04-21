from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
import numpy as np
import pandas as pd
from pathlib import Path
import pyBigWig
from typing import Iterator, Optional, Tuple
from tqdm import tqdm
from pandas.core.groupby import DataFrameGroupBy

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

    def load_coverage(self, norm_covg_dirs: list, batch_id: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
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

    def by_gene(self) -> Iterator[Tuple[str, pd.DataFrame]]:
        """Iterate over coverage data for each gene
        
        Yields:
            Tuple of gene ID and coverage DataFrame for the gene
        """
        for gene_id, df in self.coverage.groupby('gene_id'):
            yield str(gene_id), df

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
    # Use rsplit to split on the rightmost 2 underscores in case the gene ID
    # contains underscores.
    bins.index = bins.index.str.rsplit("_", n=2, expand=True)
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

def bin_covg_from_bigwigs(bigwig_manifest: pd.DataFrame, bins: pd.DataFrame, median_coverage: Optional[float] = None) -> pd.DataFrame:
    """Load binned coverage data from bigWig files and scale by sequencing depth
    
    Args:
        bigwig_manifest: DataFrame containing bigWig manifest. Must have columns
          sample and path.
        bins: DataFrame containing bin information. Rows are bins, indexed by
          gene_id and pos, and columns include chrom, chrom_start, and
          chrom_end.
        median_coverage: Median of sumData across all samples. If provided,
          coverage values will be scaled by sumData/median_coverage to normalize
          for sequencing depth. If None, no scaling is applied.

    Returns:
        The coverage DataFrame. Rows are bins and columns are sample IDs.

    Raises:
        ValueError: If any chromosome in the bins is not found in any bigWig file
    """
    covg = np.zeros((bins.shape[0], len(bigwig_manifest)))
    for i, (_, row) in enumerate(tqdm(bigwig_manifest.iterrows(), 
                                     total=len(bigwig_manifest), 
                                     desc="Processing bigWig files")):
        with pyBigWig.open(str(row['path'])) as bw:
            # Check that all required chromosomes are in the bigWig file
            chroms = bw.chroms()
            unique_chroms = set(bins['chrom'].unique())
            missing_chroms = [chr for chr in unique_chroms if str(chr) not in chroms]
            if missing_chroms:
                available_chroms = list(chroms.keys())
                raise ValueError(
                    f"Chromosomes {missing_chroms} not found in bigWig file {row['path']}.\n"
                    f"Available chromosomes: {available_chroms}"
                )

            # Get scaling factor based on bigWig header
            scaling_factor = 1.0
            if median_coverage is not None:
                total_coverage = bw.header()['sumData']
                if total_coverage > 0:
                    scaling_factor = total_coverage / median_coverage
                    
            for j, (gene_id, pos) in enumerate(bins.index):
                chrom = str(bins.loc[(gene_id, pos), 'chrom'])
                start = bins.loc[(gene_id, pos), 'chrom_start']
                end = bins.loc[(gene_id, pos), 'chrom_end']
                try:
                    bin_coverage = bw.stats(chrom, start, end, type='mean', exact=True)[0]
                    if bin_coverage is not None:
                        covg[j, i] = bin_coverage / scaling_factor
                except RuntimeError as e:
                    print(f"RuntimeError for {gene_id}, {chrom}:{start}-{end} in file {row['path']}: {e}", flush=True)
                    raise e
    samples = bigwig_manifest['sample'].tolist()
    df = pd.DataFrame(covg, index=bins.index, columns=samples)
    return df

def base_covg_from_bigwigs(bigwig_paths: list[Path], seqname: str, start: int, end: int, median_coverage: Optional[float] = None) -> np.ndarray:
    """Load base-level coverage data for one region from bigWig files

    Args:
        bigwig_paths: List of paths to bigWig files to load
        seqname: Chromosome name
        start: Start position (0-based)
        end: End position (0-based)
        median_coverage: Median of sumData across all samples. If provided,
          coverage values will be scaled by sumData/median_coverage to normalize
          for sequencing depth. If None, no scaling is applied.

    Returns:
        Array of shape (end - start, len(bigwig_paths)) with coverage data

    Raises:
        ValueError: If seqname is not found in any of the bigWig files
    """
    covg = np.zeros((end - start, len(bigwig_paths)))
    for i, path in enumerate(bigwig_paths):
        with pyBigWig.open(str(path)) as bw:
            chroms = bw.chroms()
            if str(seqname) not in chroms:
                available_chroms = list(chroms.keys())
                raise ValueError(
                    f"Chromosome '{seqname}' not found in bigWig file {path}.\n"
                    f"Available chromosomes: {available_chroms}"
                )
            
            # Get scaling factor based on bigWig header
            scaling_factor = 1.0
            if median_coverage is not None:
                total_coverage = bw.header()['sumData']
                if total_coverage > 0:
                    scaling_factor = total_coverage / median_coverage
                    
            try: # Print interval, then throw error
                values = bw.values(str(seqname), start, end)
                covg[:, i] = np.array(values) / scaling_factor if values is not None else 0
            except RuntimeError as e:
                print(f"RuntimeError for {seqname}:{start}-{end} in file {path}: {e}", flush=True)
                raise e
    return covg

def validate_chromosomes(bigwig_path: Path, required_chroms: set[str]) -> None:
    """Validate that a bigWig file contains all required chromosomes

    One issue this checks for is mismatched chromosome formats between the
    bigWig file and the annotations used for binning.

    Args:
        bigwig_path: Path to a bigWig file
        required_chroms: Set of chromosome names that must be present

    Raises:
        ValueError: If any required chromosome is not found in the bigWig file
    """
    with pyBigWig.open(str(bigwig_path)) as bw:
        chroms = bw.chroms()
        missing_chroms = [chr for chr in required_chroms if str(chr) not in chroms]
        if missing_chroms:
            available_chroms = list(chroms.keys())
            raise ValueError(
                f"Chromosomes {missing_chroms} not found in bigWig file {bigwig_path}.\n"
                f"Available chromosomes: {available_chroms}"
            )

def load_phenotypes(pheno_file: Path, samples: list) -> pd.DataFrame:
    """Load phenotype table in BED format
    
    Gene IDs are parsed from the 4th column from the start up to the first
    double underscore. Chromosome, start, and end columns are ignored.

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

def load_and_prepare_phenos(pheno_files: list, samples: list, n_pcs: int) -> DataFrameGroupBy:
    """Load explicit phenotype tables and prepare for regression

    Load phenotype tables, remove principal components, and group by gene ID.
    Zero-variance phenotypes are filtered out before grouping by gene.
    Phenotype IDs must begin with the gene ID followed by two underscores. A
    period and subsequent characters in a gene ID are assumed to be version
    information and are removed, e.g. `ENSG00000000457.14__ENST00000367771.11`
    will be matched with latent phenotypes for gene ENSG00000000457.

    Args:
        pheno_files: List of paths to phenotype tables.
        samples: List of sample IDs to filter phenotypes by.
        n_pcs: Number of principal components to remove.

    Returns:
        The prepared phenotype DataFrame grouped by gene ID. Rows are phenotypes
        and columns are samples.
    """
    print(f'Loading phenotypes from {len(pheno_files)} file{"" if len(pheno_files) == 1 else "s"} to regress out...', flush=True)
    phenos = [load_phenotypes(f, samples).reset_index() for f in pheno_files]
    phenos = pd.concat(phenos, axis=0)
    # Parse gene IDs from phenotype IDs
    pheno_genes = phenos['phenotype_id'].str.split('__').str[0]
    # Remove version information from gene IDs
    pheno_genes = pheno_genes.str.split('.').str[0]
    # Phenotype IDs can't be used as index due to duplicates across files
    phenos = phenos.drop('phenotype_id', axis=1)
    # Remove covariates using PCA
    if n_pcs > 0:
        phenos = remove_pcs_from_phenos(phenos, n_pcs)
    # Filter out zero-variance features (rows)
    phenos = phenos.loc[phenos.var(axis=1) > 0]
    phenos = phenos.groupby(pheno_genes)
    return phenos

def regress_out_phenos_single(y: pd.Series, x: np.ndarray) -> pd.Series:
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
    x_const = sm.add_constant(x, has_constant='add')
    model_all = sm.OLS(y, x_const).fit()
    assert model_all.pvalues.index[0] == 'const'
    x_sig = x[:, model_all.pvalues[1:] < 0.01] # exclude intercept
    if x_sig.shape[1] == 0:
        return y
    model = LinearRegression().fit(x_sig, y)
    resid = y - model.predict(x_sig)
    return pd.Series(resid, index=y.index)

def regress_out_phenos_gene(df: pd.DataFrame, phenos: DataFrameGroupBy, gene_id: str) -> pd.DataFrame:
    """Regress out phenotypes per gene
    
    Args:
        df: The input DataFrame containing normalized coverage per bin per
          sample for one gene.
        phenos: The input DataFrameGroupBy containing phenotypes to regress out.
          Rows are phenotypes and columns are samples.
        gene_id: The gene ID to process.

    Returns:
        The DataFrame of residuals after regressing out phenotypes.
    """
    if gene_id in phenos.groups:
        x = phenos.get_group(gene_id).to_numpy().T
        df = df.apply(regress_out_phenos_single, axis=1, args=(x,)) # axis=1 means by row, i.e. index of Series is column names
    return df

def remove_pcs_from_phenos(phenos: pd.DataFrame, n_pcs: int) -> pd.DataFrame:
    """Remove principal components from phenotype data

    This corrects for technical covariates and other confounding factors
    across phenotypes, and keeping more gene-specific variation, before
    regressing the phenotypes out of the coverage data.
    
    Args:
        phenos: DataFrame of phenotypes (rows) x samples (columns)
        n_pcs: Number of principal components to remove. If 0, no PCs are removed.
        
    Returns:
        DataFrame of corrected phenotype values with same dimensions as input
    """
    print(f'Removing {n_pcs} principal components from phenotype data...', flush=True)
        
    # If we have fewer samples than PCs, reduce the number of PCs
    if phenos.shape[1] < n_pcs:
        n_pcs = phenos.shape[1] - 1
        print(f"Warning: Reducing n_pcs to {n_pcs} because we have only {phenos.shape[1]} samples", flush=True)
        if n_pcs <= 0:
            print("Cannot perform PCA with too few samples. Returning original data.", flush=True)
            return phenos

    # Perform PCA on samples
    pca = PCA(n_components=n_pcs)
    pc_matrix = pca.fit_transform(phenos.T)  # samples × n_pcs
    
    # For each phenotype, regress out the PCs
    corrected_phenos = phenos.copy()
    
    # This approach regresses each phenotype separately against the PCs
    for i in range(phenos.shape[0]):
        # Get phenotype values across all samples (should be 1D array with length = number of samples)
        y = phenos.iloc[i].values
        
        # Safety check: ensure dimensions match
        if len(y) != pc_matrix.shape[0]:
            print(f"Warning: Dimension mismatch for phenotype {i}. " 
                  f"Phenotype values: {len(y)}, PC matrix rows: {pc_matrix.shape[0]}")
            
        # Fit linear model: phenotype ~ PCs
        model = LinearRegression().fit(pc_matrix, y)
        
        # Calculate residuals (phenotype values after removing PC effect)
        pc_effects = model.predict(pc_matrix)
        corrected_phenos.iloc[i] = y - pc_effects
    
    print(f'Removed {n_pcs} PCs explaining {pca.explained_variance_ratio_.sum():.2%} of variance', flush=True)
    
    return corrected_phenos

def regress_out_phenos(covg: pd.DataFrame, phenos: DataFrameGroupBy) -> pd.DataFrame:
    """Regress out phenotypes from coverage data
    
    Args:
        covg: The input DataFrame containing normalized coverage per bin per
          sample for multiple genes.
        phenos: DataFrameGroupBy of phenotypes (rows) x samples (columns),
          grouped by gene ID.
    Returns:
        The DataFrame of residuals after regressing out phenotypes.
    """
    prev_index = covg.index.copy()
    
    # Track genes with explicit phenotypes
    genes_with_phenos = 0
    total_genes = 0
    
    def regress_with_tracking(c):
        nonlocal genes_with_phenos, total_genes
        total_genes += 1
        if c.name in phenos.groups:
            genes_with_phenos += 1
        return regress_out_phenos_gene(c, phenos, c.name)
    
    covg = covg.groupby("gene_id", group_keys=False).apply(regress_with_tracking)
    assert covg.index.equals(prev_index)
    
    print(f"Regression statistics: {genes_with_phenos} out of {total_genes} genes ({genes_with_phenos/total_genes:.1%}) had explicit phenotypes to regress out", flush=True)
    
    return covg

def prepare_coverage(
        bigwig_manifest: pd.DataFrame,
        bins_dir: Path,
        batches: list,
        pheno_files: list,
        outdir: Path,
        median_coverage: Optional[float] = None,
        n_pcs: int = 4,
) -> None:
    """Assemble per-sample coverage, split into batches, and save to numpy binary files
    
    Args:
        bigwig_manifest: DataFrame containing bigWig manifest. Must have columns
          sample and path.
        bins_dir: Directory of per-batch BED files containing bin regions. Must
          have start, end and bin ID in 2nd, 3rd, and 4th columns. Rows must
          correspond to rows of corresponding input coverage file.
        batches: List of batch IDs to process. Batch IDs are integers starting
          from 0.
        pheno_files: List of paths to phenotype tables to regress out of the
          coverage data per gene prior to model input.
        outdir: Directory where per-batch numpy binary files with normalized
          coverage will be written.
        median_coverage: Median of sumData across all samples. If provided,
          coverage values will be scaled by sumData/median_coverage to normalize
          for sequencing depth. If None, no scaling is applied.
        n_pcs: Number of principal components to remove from phenotype data to
          account for technical covariates. If 0, no PCs are removed.
    """
    if len(pheno_files) > 0:
        samples = bigwig_manifest['sample'].tolist()
        phenos = load_and_prepare_phenos(pheno_files, samples, n_pcs)

    for batch in batches:
        print(f'=== Preparing coverage for batch {batch} ===', flush=True)
        binfile = bins_dir / f'batch_{batch}.bed.gz'
        bins = load_bins(binfile)
        
        covg = bin_covg_from_bigwigs(bigwig_manifest, bins, median_coverage)
        # Sort by gene and position along gene instead of genomic coordinate
        covg = covg.sort_index()
        assert not covg.isna().any().any()
        # Enable log-transformation and reduce impact of noise
        covg += 8
        # Variance stabilization
        covg = np.log2(covg)
        
        if len(pheno_files) > 0:
            covg = regress_out_phenos(covg, phenos)
        
        outdir.mkdir(exist_ok=True)
        mat = covg.to_numpy()
        np.save(outdir / f'batch_{batch}.npy', mat)
        bins = bins.loc[covg.index, :]
        bins.to_csv(outdir / f'batch_{batch}.bins.tsv.gz', sep='\t')
        
        # If samples.txt doesn't exist, write it:
        if not (outdir / 'samples.txt').exists():
            with open(outdir / 'samples.txt', 'w') as f:
                f.write('\n'.join(covg.columns) + '\n')
        
        print(f'Coverage saved in {outdir}', flush=True)

