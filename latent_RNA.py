"""Extract latent transcriptomic phenotypes from RNA-seq coverage data"""

import argparse
import math
from pathlib import Path
import pickle
from typing import Iterator
import numpy as np
import pandas as pd
import pyBigWig
from skfda.preprocessing.dim_reduction import FPCA
from skfda.representation.basis import BSplineBasis
from skfda.representation.grid import FDataGrid
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA
import statsmodels.api as sm
from tqdm import tqdm

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

class Model:
    def __init__(self, fpca: bool, fpca_x_values: str, fpca_basis: str):
        """Initialize a model object

        Args:
            fpca: Whether to fit a functional PCA model instead of regular PCA.
            fpca_x_values: Whether to use bin numbers or genomic positions as
              x-values for functional PCA. Options are 'bin' and 'pos'.
            fpca_basis: Basis function to use for functional PCA. Options are
              'discrete' for discretized FPCA directly on the data, and 'spline'
              for a 4th-order B-spline basis.
        """
        self.model = None
        self.features = None
        self.fpca = fpca
        self.fpca_x_values = fpca_x_values
        self.fpca_basis = fpca_basis

    def fit(self, df: pd.DataFrame, var_expl: float, n_pcs_max: int):
        """Fit functional PCA model to normalized coverage for one gene
        
        Saves enough PCs to explain var_expl variance or n_pcs_max, whichever
        is fewer.

        Args:
            df: The input DataFrame containing normalized coverage per bin per
              sample.
            var_expl: Maximum variance explained by the PCs kept per gene. Pass
              0 or 1 for no variance explained cutoff.
            n_pcs_max: Maximum number of PCs to keep per gene. Pass 0 for no
              cutoff.
        """
        if self.fpca:
            n_bins_with_var = (df.std(axis=1) > 0).sum()
            x = df.values.T
            if self.fpca_x_values == 'bin':
                fd = FDataGrid(x, grid_points=np.arange(x.shape[1]))
            else:
                fd = FDataGrid(x, grid_points=df.index.get_level_values('pos'))
            # FPCA currently defaults to n_components=3, so we need to set it explicitly:
            # (There's no option for variance explained cutoff, so we'll subset the PCs later.)
            if n_pcs_max == 0:
                n_pcs_max = x.shape[0]
            n_comp = min(n_pcs_max, x.shape[0], n_bins_with_var)
            if n_comp == 0:
                return
            if self.fpca_basis == 'spline':
                self.SPLINE_ORDER = 4
                if n_comp < self.SPLINE_ORDER:
                    return
                fd = fd.to_basis(BSplineBasis(n_basis=n_comp, order=self.SPLINE_ORDER))
            # Must pass weights to avoid 'Matrix is not positive definite' error:
            model = FPCA(n_components=n_comp, _weights=np.ones(x.shape[1]))
            try:
                model.fit(fd)
            except Exception as e:
                print(f"Error occurred during model fitting: {e}", flush=True)
                return
            if var_expl not in {0, 1}:
                subset_pcs_fpca(model, var_expl, self.fpca_basis)
        else:
            df = df.loc[df.std(axis=1) > 0, :]
            if df.shape[0] == 0:
                return
            x = df.values.T
            # Center and scale the data:
            xmean = x.mean(axis=0)
            xstd = x.std(axis=0)
            x = (x - xmean) / xstd
            model = PCA(n_components=var_expl) if var_expl not in {0, 1} else PCA()
            try:
                model.fit(x)
            except np.linalg.LinAlgError:
                print("SVD did not converge, retrying with added noise", flush=True)
                x += np.random.normal(0, 1e-6, size=x.shape)
                model.fit(x)
            subset_pcs_pca(model, n_pcs_max)
            self.features = pd.DataFrame(index=df.index) # Save feature stats for transforming new data
            self.features['mean'] = xmean
            self.features['std'] = xstd
        self.model = model

    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply functional PCA model to normalized coverage for one gene
        
        Args:
            df: The input DataFrame containing normalized coverage per bin per
              sample.

        Returns:
            The DataFrame of transformed data, with gene ID and PC number as
            index levels and sample IDs as columns.
        """
        # Filter df to include same features as the model:
        if self.fpca:
            x = df.values.T
            if self.fpca_x_values == 'bin':
                fd = FDataGrid(x, grid_points=np.arange(x.shape[1]))
            else:
                fd = FDataGrid(x, grid_points=df.index.get_level_values('pos'))
            if self.fpca_basis == 'spline':
                n_basis = self.model.components_.n_basis
                fd = fd.to_basis(BSplineBasis(n_basis=n_basis, order=self.SPLINE_ORDER))
            mat = self.model.transform(fd).T
        else:
            df = df.loc[self.features.index, :]
            x = df.values.T
            # Center and scale the data using stats from the model-fitting data:
            x = (x - self.features['mean'].values) / self.features['std'].values
            mat = self.model.transform(x).T
        pc_names = [f'PC{i + 1}' for i in range(mat.shape[0])]
        out = pd.DataFrame(mat, index=pc_names, columns=df.columns)
        out.index.set_names('PC', inplace=True)
        # Add gene_id to index:
        out.reset_index(inplace=True)
        out.insert(0, 'gene_id', df.index.get_level_values('gene_id')[0])
        out.set_index(['gene_id', 'PC'], inplace=True)
        return out

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

def covg_from_bigwigs(bigwig_paths_file: Path, bins: pd.DataFrame) -> pd.DataFrame:
    """Load coverage data from bigWig files
    
    Args:
        bigwig_paths_file: Path to a file containing paths to bigWig files.
          Sample IDs will be extracted from the basenames of the files.
        bins: DataFrame containing bin information. Rows are bins, indexed by
          gene_id and pos, and columns include chrom, chrom_start, and
          chrom_end.

    Returns:
        The coverage DataFrame. Rows are bins and columns are sample IDs.
    """
    with open(bigwig_paths_file) as f:
        bigwig_paths = [Path(l.strip()) for l in f]
    covg = np.zeros((bins.shape[0], len(bigwig_paths)))
    for i, bw in enumerate(bigwig_paths):
        with pyBigWig.open(str(bw)) as f:
            for j, (gene_id, pos) in enumerate(bins.index):
                chrom = bins.loc[(gene_id, pos), 'chrom']
                start = bins.loc[(gene_id, pos), 'chrom_start']
                end = bins.loc[(gene_id, pos), 'chrom_end']
                try:
                    covg[j, i] = f.stats(chrom, start, end, type='mean', exact=True)[0]
                except RuntimeError as e:
                    print(f"RuntimeError for {gene_id}, {chrom}:{start}-{end} in file {bw.stem}: {e}", flush=True)
    samples = [bw.stem for bw in bigwig_paths]
    df = pd.DataFrame(covg, index=bins.index, columns=samples)
    return df

def normalize_coverage(df: pd.DataFrame, bins: pd.DataFrame, pseudocount: float = 8) -> pd.DataFrame:
    """Normalize coverage for each gene in one batch.

    1. A pseudocount is added to the mean coverage to allow log-transformation
       and reduce the impact of noise in low-coverage bins.
    2. Coverage is normalized to fraction of gene's coverage for each bin, to
       control for total gene expression and sequencing depth.
    3. Coverage is normalized to mean fraction of gene's coverage per bp for
       each bin.

    Args:
        df: The input DataFrame containing mean coverage per bin per sample.
          Rows are bins and columns are samples. Multi-indexed by gene_id and
          pos.
        bins: The input DataFrame containing bin information.
        pseudocount: Pseudocount to add to coverage values (mean coverage per
          bp for each bin) before normalization.

    Returns:
        The normalized coverage DataFrame.

    Raises:
        AssertionError: If the input DataFrame contains any NaN values.
    """
    assert df.index.equals(bins.index)
    # Add pseudocount to mean coverage per bp for each bin
    df += pseudocount
    # Convert to total coverage per bin
    lengths = bins['chrom_end'] - bins['chrom_start']
    df = df.mul(lengths, axis=0)
    # Normalize to fraction of gene's coverage for each bin
    df = df.groupby('gene_id', group_keys=False).apply(lambda x: x.div(x.sum(axis=0), axis=1))
    assert not df.isna().any().any()
    # Normalize to mean fraction of gene's coverage per bp for each bin
    df = df.div(lengths, axis=0)
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

def prepare(bigwig_paths_file: Path, bins_dir: Path, batches: list, pheno_files: list, outdir: Path):
    """Assemble per-sample coverage, split into batches, and save to numpy binary files
    
    Args:
        bigwig_paths_file: Path to a file containing paths to bigWig files.
          Sample IDs will be extracted from the basenames of the files.
        bins_dir: Directory of per-batch BED files containing bin regions. Must
          have start, end and bin ID in 2nd, 3rd, and 4th columns. Rows must
          correspond to rows of corresponding input coverage file.
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
        covg = covg_from_bigwigs(bigwig_paths_file, bins)
        # Sort by gene and position along gene instead of genomic coordinate
        covg = covg.sort_index()
        bins = bins.loc[covg.index, :]
        covg = normalize_coverage(covg, bins)
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

def subset_pcs_fpca(fpca: FPCA, var_expl: float, fpca_basis: str):
    """Subset the FPCA model using variance explained cutoff
    
    Args:
        fpca: The input FPCA model.
        var_expl: Maximum variance explained by the PCs kept per gene.
        fpca_basis: Basis function to use for functional PCA. Options are
          'discrete' for discretized FPCA directly on the data, and 'spline'
          for a 4th-order B-spline basis.
    """
    # Find the number of PCs needed to explain var_expl variance:
    n_pcs_max = np.argmax(np.cumsum(fpca.explained_variance_ratio_) >= var_expl) + 1
    if fpca_basis == 'discrete':
        if fpca.components_.shape[0] <= n_pcs_max:
            return
        fpca.components_.data_matrix = fpca.components_.data_matrix[:n_pcs_max, :]
    else:
        if fpca.components_.coefficients.shape[0] <= n_pcs_max:
            return
        fpca.components_.coefficients = fpca.components_.coefficients[:n_pcs_max, :]
    fpca.explained_variance_ = fpca.explained_variance_[:n_pcs_max]
    fpca.explained_variance_ratio_ = fpca.explained_variance_ratio_[:n_pcs_max]
    fpca.singular_values_ = fpca.singular_values_[:n_pcs_max]

def subset_pcs_pca(pca: PCA, n_pcs_max: int):
    """Subset the PCA model to n_pcs_max PCs to reduce file size
    
    Args:
        pca: The input PCA model.
        n_pcs_max: Maximum number of PCs to keep per gene. Pass 0 for no cutoff.
    """
    if n_pcs_max != 0 and n_pcs_max < pca.n_components_:
        pca.components_ = pca.components_[:n_pcs_max, :]
        pca.n_components_ = n_pcs_max
        pca.explained_variance_ = pca.explained_variance_[:n_pcs_max]
        pca.explained_variance_ratio_ = pca.explained_variance_ratio_[:n_pcs_max]
        pca.singular_values_ = pca.singular_values_[:n_pcs_max]
        pca.noise_variance_ = None # Could compute if necessary, but omitting otherwise

def fit_batch(norm_covg_dirs: list, batch: int, var_expl_max: float, n_pcs_max: int, fpca: bool = False, fpca_x_values: str = 'bin', fpca_basis: str = 'discrete') -> dict:
    """Fit functional PCA models for all genes in a batch
    
    Args:
        norm_covg_dirs: List of directories containing normalized coverage data.
        batch: ID of the batch to fit.
        var_expl_max: Maximum variance explained by the PCs kept per gene. Pass
          0 or 1 for no variance explained cutoff.
        n_pcs_max: Maximum number of PCs to keep per gene. Pass 0 for no cutoff.
        fpca: Whether to fit a functional PCA model instead of regular PCA.
        fpca_x_values: Whether to use bin numbers or genomic positions as
          x-values for functional PCA. Options are 'bin' and 'pos'.
        fpca_basis: Basis function to use for functional PCA. Options are
          'discrete' for discretized FPCA directly on the data, and 'spline' for
          a 4th-order B-spline basis.

    Returns:
        Dictionary of fitted models, with keys 'var_expl_max', 'n_pcs_max', and
        'models'. The 'models' key contains a dictionary with gene IDs as keys
        and Model objects as values.
    """
    covg = CoverageData(norm_covg_dirs, batch)
    models = {'var_expl_max': var_expl_max, 'n_pcs_max': n_pcs_max, 'models': {}}
    for gene_id, x in tqdm(covg.by_gene(), total=len(covg.genes), desc="Fitting models"):
        model = Model(fpca, fpca_x_values, fpca_basis)
        model.fit(x, var_expl_max, n_pcs_max)
        if model.model is not None:
            models['models'][gene_id] = model
    return models

def fit(norm_covg_dirs: list, batches: list, var_expl_max: float, n_pcs_max: int, output_dir: Path, fpca: bool = False, fpca_x_values: str = 'bin', fpca_basis: str = 'discrete'):
    """Fit functional PCA models for all genes in one or more batches
    
    Args:
        norm_covg_dirs: List of directories containing normalized coverage data.
        batches: List of batch IDs to fit.
        var_expl_max: Maximum variance explained by the PCs kept per gene. Pass
          0 or 1 for no variance explained cutoff.
        n_pcs_max: Maximum number of PCs to keep per gene. Pass 0 for no cutoff.
        output_dir: Directory in which to save model pickle files.
        fpca: Whether to fit a functional PCA model instead of regular PCA.
        fpca_x_values: Whether to use bin numbers or genomic positions as
          x-values for functional PCA. Options are 'bin' and 'pos'.
        fpca_basis: Basis function to use for functional PCA. Options are
          'discrete' for discretized FPCA directly on the data, and 'spline' for
          a 4th-order B-spline basis.
    """
    output_dir.mkdir(exist_ok=True)
    for batch in batches:
        print(f'=== Fitting models for batch {batch} ===', flush=True)
        models = fit_batch(norm_covg_dirs, batch, var_expl_max, n_pcs_max, fpca, fpca_x_values, fpca_basis)
        outfile = output_dir / f'models_batch_{batch}.pickle'
        with open(outfile, 'wb') as f:
            pickle.dump(models, f)
        print(f'Models saved to {outfile}', flush=True)

def transform_batch(norm_covg_dir: Path, batch: int, models_dir: dict) -> Iterator[pd.DataFrame]:
    """Apply functional PCA transformation to all genes in a batch
    
    Args:
        norm_covg_dir: Directory of per-batch numpy binary files with normalized
          coverage.
        batch: ID of the batch to transform.
        models_dir: Directory of saved models to load and use for transformation.

    Yields:
        DataFrame of transformed data for one gene, with gene ID and PC number
        as index levels and sample IDs as columns.
    """
    covg = CoverageData([norm_covg_dir], batch)
    print('  Loading models...', flush=True)
    models_file = models_dir / f'models_batch_{batch}.pickle'
    with open(models_file, 'rb') as f:
        models = pickle.load(f)
    for gene_id, x in tqdm(covg.by_gene(), total=len(covg.genes), desc="  Generating phenotypes per gene"):
        if gene_id in models['models']:
            yield models['models'][gene_id].transform(x)

def transform(norm_covg_dir: Path, models_dir: dict, n_batches: int) -> pd.DataFrame:
    """Apply functional PCA transformation to all genes
    
    Args:
        norm_covg_dir: Directory of per-batch numpy binary files with normalized
          coverage.
        models_dir: Directory of saved models to load and use for transformation.
        n_batches: Number of batches in the data. Latent phenotypes from all
          batches will be computed and concatenated.
    """
    out = []
    for batch in range(n_batches):
        print(f'Transforming batch {batch}', flush=True)
        out.extend(transform_batch(norm_covg_dir, batch, models_dir))
    out = pd.concat(out)
    out = out.reset_index()
    return out

def create_parser():
    parser = argparse.ArgumentParser(description='Fit and/or apply models on feature bin coverage data')
    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand', required=True, help='Choose a subcommand')

    # Subparser for 'prepare'
    parser_prepare = subparsers.add_parser('prepare', help='Prepare RNA-seq coverage for a batch of genes')
    parser_prepare.add_argument('-i', '--bigwig-paths-file', type=Path, metavar='FILE', required=True, help='File containing list of paths to per-sample bigWig files. Basenames of files will be used as sample IDs.')
    parser_prepare.add_argument('--bins-dir', type=Path, metavar='DIR', required=True, help='Directory of per-batch BED files containing bin regions.')
    parser_prepare_batch = parser_prepare.add_mutually_exclusive_group(required=True)
    parser_prepare_batch.add_argument('-b', '--batch', type=int, metavar='N', help='Batch ID to process. Batch IDs are integers starting from 0.')
    parser_prepare_batch.add_argument('--n-batches', type=int, metavar='N', help='To load and fit all batches in sequence, provide number of batches instead of a specific batch.')
    parser_prepare_phenos = parser_prepare.add_mutually_exclusive_group(required=False)
    parser_prepare_phenos.add_argument('-p', '--pheno-paths', nargs='+', type=Path, metavar='FILE', help='One or more paths to phenotype tables to regress out of the coverage data per gene prior to model input. Files should be in bed format, i.e. input format for tensorqtl. Gene IDs are parsed from the 4th column from the start up to the first non-alphanumeric character.')
    parser_prepare_phenos.add_argument('--pheno-paths-file', type=Path, metavar='FILE', help='File containing list of paths to phenotype tables.')
    parser_prepare.add_argument('-o', '--output-dir', type=Path, metavar='DIR', required=True, help='Directory where per-batch numpy binary files with normalized coverage will be written.')

    # Subparser for 'fit'
    parser_fit = subparsers.add_parser('fit', help='Fit FPCA or PCA models to coverage data')
    parser_fit_input = parser_fit.add_mutually_exclusive_group(required=True)
    parser_fit_input.add_argument('-d', '--norm-covg-dir', nargs='+', type=Path, metavar='DIR', help='Directory of per-batch numpy binary files with normalized coverage. Specify multiple directories to load data from all of them and fit models using the combined dataset. All datasets must have been generated using the same per-batch gene bins files.')
    parser_fit_input.add_argument('--norm-covg-dir-file', type=Path, metavar='FILE', help='File containing list of directories of per-batch numpy binary files. Use this instead of -d/--norm-covg-dir in case of many directories.')
    parser_fit_batch = parser_fit.add_mutually_exclusive_group(required=True)
    parser_fit_batch.add_argument('-b', '--batch', type=int, metavar='N', help='Batch ID to load and fit. Batch IDs are integers starting from 0.')
    parser_fit_batch.add_argument('--n-batches', type=int, metavar='N', help='To load and fit all batches in sequence, provide number of batches instead of a specific batch.')
    parser_fit.add_argument('-m', '--models-dir', type=Path, metavar='DIR', required=True, help='Directory in which to save model pickle files.')
    parser_fit.add_argument('-v', '--var-expl-max', type=float, default=0.8, metavar='FLOAT', help='Max variance explained by the PCs kept per gene. Pass 0 or 1 for no variance explained cutoff.  (default: %(default)s)')
    parser_fit.add_argument('-n', '--n-pcs-max', type=int, default=16, metavar='N', help='Max number of PCs to keep per gene. Pass 0 for no cutoff. (default: %(default)s)')
    parser_fit.add_argument('--use-fpca', action='store_true', help='Use functional PCA instead of regular PCA.')
    parser_fit.add_argument('--fpca-x-values', type=str, choices=['bin', 'pos'], default='bin', metavar='STRING', help='Whether to use bin numbers or genomic positions as x-values for functional PCA. (default: %(default)s)')
    parser_fit.add_argument('--fpca-basis', type=str, choices=['discrete', 'spline'], default='discrete', metavar='STRING', help='Basis function to use for functional PCA. `discrete` will run discretized FPCA directly on the data, `spline` will use a 4th-order B-spline basis. (default: %(default)s)')

    # Subparser for 'transform'
    parser_transform = subparsers.add_parser('transform', help='Apply fitted models to coverage data')
    parser_transform.add_argument('-d', '--norm-covg-dir', type=Path, metavar='DIR', required=True, help='Directory of per-batch numpy binary files with normalized coverage.')
    parser_transform.add_argument('-m', '--models-dir', type=Path, metavar='DIR', required=True, help='Directory of saved models (*.pickle) to load and use for transformation.')
    parser_transform.add_argument('--n-batches', type=int, metavar='N', required=True, help='Number of batches in the data. Latent phenotypes from all batches will be computed and concatenated.')
    parser_transform.add_argument('-o', '--output', type=Path, metavar='FILE', required=True, help='Output file (TSV).')

    return parser

if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()
    if args.subcommand == 'prepare':
        if args.pheno_paths is not None:
            pheno_paths = args.pheno_paths
        elif args.pheno_paths_file is not None:
            with open(args.pheno_paths_file) as f:
                pheno_paths = [Path(l.strip()) for l in f]
        else:
            pheno_paths = []
        batches = [args.batch] if args.batch is not None else list(range(args.n_batches))
        prepare(args.bigwig_paths_file, args.bins_dir, batches, pheno_paths, args.output_dir)
    elif args.subcommand == 'fit':
        if args.norm_covg_dir is not None:
            norm_covg_dirs = args.norm_covg_dir
        else:
            with open(args.norm_covg_dir_file) as f:
                norm_covg_dirs = [Path(l.strip()) for l in f]
        batches = [args.batch] if args.batch is not None else list(range(args.n_batches))
        fit(norm_covg_dirs, batches, args.var_expl_max, args.n_pcs_max, args.models_dir, args.use_fpca, args.fpca_x_values, args.fpca_basis)
    elif args.subcommand == 'transform':
        print('=== Generating latent phenotypes ===', flush=True)
        output = transform(args.norm_covg_dir, args.models_dir, args.n_batches)
        output.to_csv(args.output, sep='\t', index=False, float_format='%g')
        print(f'Latent phenotypes saved to {args.output}', flush=True)
