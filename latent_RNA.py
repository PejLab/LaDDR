"""Extract latent transcriptomic phenotypes from RNA-seq coverage data"""

import argparse
import math
from pathlib import Path
import pickle
from typing import Iterator
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from skfda.preprocessing.dim_reduction import FPCA
from skfda.representation.basis import BSplineBasis
from skfda.representation.grid import FDataGrid
from sklearn.linear_model import LinearRegression
from tqdm import tqdm

class CoverageData:
    def __init__(self, batch_covg_dirs: list, batch_id: int):
        self.samples = []
        for d in batch_covg_dirs:
            with open(d / 'samples.txt') as f:
                self.samples.extend([l.strip() for l in f])
        if len(batch_covg_dirs) > 1:
            # Make sample IDs unique for concatenation, and ID won't be saved in the models anyway.
            self.samples = [f'S{i + 1}' for i in range(len(self.samples))]
        self.coverage, self.regions = self.load_coverage(batch_covg_dirs, batch_id)
        self.genes = list(self.regions.groupby('gene_id').groups.keys())

    def load_coverage(self, batch_covg_dirs: list, batch_id: int) -> tuple:
        mats = [np.load(d / f'covg_{batch_id}.npy') for d in batch_covg_dirs]
        mat = np.concatenate(mats, axis=1)
        region_file = batch_covg_dirs[0] / f'covg_{batch_id}.regions.tsv.gz'
        regions = pd.read_csv(region_file, sep='\t', index_col=[0, 1])
        assert mat.shape[0] == regions.shape[0]
        assert mat.shape[1] == len(self.samples)
        df = pd.DataFrame(mat, index=regions.index, columns=self.samples)
        return df, regions

    def by_gene(self) -> Iterator[tuple]:
        for gene_id, df in self.coverage.groupby('gene_id'):
            yield gene_id, df

class Model:
    def __init__(self, fpca: bool, fpca_x_values: str, fpca_basis: str):
        self.model = None
        self.features = None
        self.fpca = fpca
        self.fpca_x_values = fpca_x_values
        self.fpca_basis = fpca_basis

    def fit(self, df: pd.DataFrame, var_expl: float, n_pcs_max: int):
        """Fit functional PCA model to normalized coverage for one gene
        
        Saves enough PCs to explain var_expl variance or n_pcs_max, whichever is smaller.
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
        """Apply functional PCA model to normalized coverage for one gene"""
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

def load_regions(regionfile: Path) -> pd.DataFrame:
    regions = pd.read_csv(regionfile, sep='\t', header=None, usecols=[3], names=['region'], index_col='region')
    regions.index = regions.index.str.split('_', n=2, expand=True)
    regions.index.set_names(['gene_id', 'start', 'end'], inplace=True)
    regions = regions.reset_index()
    regions['start'] = regions['start'].astype(int)
    regions['end'] = regions['end'].astype(int)
    regions['length'] = regions['end'] - regions['start']
    regions['pos'] = regions['start'] + regions['length'] // 2
    regions = regions.set_index(['gene_id', 'pos'])
    regions = regions.drop(['start', 'end'], axis=1)
    return regions

def load_sample_covg_np(infiles: list, nrow: int) -> np.ndarray:
    """Load per-sample coverage files and concatenate columns into one numpy array"""
    counts = np.zeros((nrow, len(infiles)), dtype=np.int32)
    for i, fname in tqdm(enumerate(infiles), total=len(infiles), desc="Loading per-sample coverage files"):
        d = pd.read_csv(fname, header=None, dtype=np.int32, names=['count'], index_col=False)
        assert d.shape[0] == nrow
        counts[:, i] = d['count']
    return counts

def load_sample_covg(infile: Path, regions: pd.DataFrame) -> pd.DataFrame:
    """Load per-sample coverage files into one dataframe"""
    infiles = pd.read_csv(infile, header=None, names=['path'])
    infiles['sample'] = [Path(p).name.split('.')[0] for p in infiles['path']]
    assert infiles['sample'].is_unique
    counts = load_sample_covg_np(infiles['path'], nrow=regions.shape[0])
    counts = pd.DataFrame(counts, index=regions.index, columns=infiles['sample'])
    # Sort by gene and position along gene instead of genomic coordinate
    counts = counts.sort_index()
    return counts

def normalize_coverage(df: pd.DataFrame, regions: pd.DataFrame) -> pd.DataFrame:
    """Normalize coverage for each gene in one batch"""
    # For log-transform, add pseudocount of 1 per bp:
    df = df.add(regions['length'].loc[df.index], axis=0)
    # Normalize coverage to ignore total expression, i.e. total coverage per sample per gene sums to 1:
    df = df.groupby('gene_id', group_keys=False).apply(lambda x: x.div(x.sum(axis=0), axis=1))
    assert not df.isna().any().any()
    # Normalize to mean read fraction per bp within each region:
    df = df.div(regions['length'].loc[df.index], axis=0)
    # Variance stabilization:
    df = np.log2(df)
    return df

def load_phenotypes(pheno_file: Path, samples: list) -> pd.DataFrame:
    """Load phenotype table in BED format
    
    Gene IDs are parsed from the 4th column from the start up to the first
    non-alphanumeric character. Chromosome, start, and end columns are ignored.
    """
    df = pd.read_csv(pheno_file, sep='\t', index_col=3, dtype={0: str, 1: str, 2: str})
    df = df.drop(df.columns[:3], axis=1)
    df.index.name = 'phenotype_id'
    # Check that all samples are present:
    missing = set(samples) - set(df.columns)
    assert len(missing) == 0, f'Phenotype table is missing samples: {missing}'
    df = df[samples]
    return df

def regress_out_phenos(df: pd.DataFrame, phenos: pd.core.groupby.DataFrameGroupBy, gene_id: str) -> pd.DataFrame:
    """Regress out phenotypes per gene"""
    if gene_id in phenos.groups:
        x = phenos.get_group(gene_id).to_numpy().T
        y = df.to_numpy().T
        model = LinearRegression().fit(x, y)
        y = y - model.predict(x)
        df = pd.DataFrame(y.T, index=df.index, columns=df.columns)
    return df

def prepare(infile: Path, regionfile: Path, pheno_files: list, batch_size: int, outdir: Path):
    """Assemble per-sample coverage, split into batches, and save to numpy binary files"""
    regions = load_regions(regionfile)
    counts = load_sample_covg(infile, regions)
    samples = counts.columns
    genes = regions.index.get_level_values('gene_id').unique().sort_values()

    n_batches = math.ceil(len(genes) / batch_size)
    gene_batch = {}
    for i, gene_id in enumerate(genes):
        gene_batch[gene_id] = i // batch_size
    regions['batch'] = [gene_batch[gene_id] for gene_id in regions.index.get_level_values('gene_id')]

    outdir.mkdir(exist_ok=True)
    with open(outdir / 'samples.txt', 'w') as f:
        f.write('\n'.join(samples) + '\n')
    with open(outdir / 'genes.txt', 'w') as f:
        f.write('\n'.join(genes) + '\n')
    with open(outdir / 'batches.txt', 'w') as f:
        f.write('\n'.join([str(i) for i in range(n_batches)]) + '\n')
    batches = counts.groupby(regions['batch'])
    regions = regions.drop('batch', axis=1)
    if len(pheno_files) > 0:
        print(f'Loading phenotypes from {len(pheno_files)} file{"" if len(pheno_files) == 1 else "s"} to regress out...', flush=True)
        phenos = [load_phenotypes(f, samples).reset_index() for f in pheno_files]
        phenos = pd.concat(phenos, axis=0)
        # Parse gene IDs from phenotype IDs:
        pheno_genes = phenos['phenotype_id'].str.extract(r'([a-zA-Z0-9]+)')[0]
        phenos = phenos.drop('phenotype_id', axis=1)
        phenos = phenos.groupby(pheno_genes)
    for batch_id, x in tqdm(batches, total=n_batches, desc="Processing and saving per-batch coverage files"):
        x = normalize_coverage(x, regions)
        if len(pheno_files) > 0:
            # Regress out phenotypes:
            prev_index = x.index.copy()
            x = x.groupby('gene_id', group_keys=False).apply(lambda c: regress_out_phenos(c, phenos, c.name))
            assert x.index.equals(prev_index)
        mat = x.to_numpy()
        np.save(outdir / f'covg_{batch_id}.npy', mat)
        reg = regions.loc[x.index, :]
        reg.to_csv(outdir / f'covg_{batch_id}.regions.tsv.gz', sep='\t')

def subset_pcs_fpca(fpca: FPCA, var_expl: float, fpca_basis: str):
    """Subset the FPCA model using variance explained cutoff"""
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
    """Subset the PCA model to n_pcs_max PCs to reduce file size"""
    if n_pcs_max != 0 and n_pcs_max < pca.n_components_:
        pca.components_ = pca.components_[:n_pcs_max, :]
        pca.n_components_ = n_pcs_max
        pca.explained_variance_ = pca.explained_variance_[:n_pcs_max]
        pca.explained_variance_ratio_ = pca.explained_variance_ratio_[:n_pcs_max]
        pca.singular_values_ = pca.singular_values_[:n_pcs_max]
        pca.noise_variance_ = None # Could compute if necessary, but omitting otherwise

def fit_batch(batch_covg_dirs: list, batch: int, var_expl_max: float, n_pcs_max: int, fpca: bool = True, fpca_x_values: str = 'bin', fpca_basis: str = 'discrete') -> dict:
    """Fit functional PCA models for all genes in a batch"""
    covg = CoverageData(batch_covg_dirs, batch)
    models = {'var_expl_max': var_expl_max, 'n_pcs_max': n_pcs_max, 'models': {}}
    for gene_id, x in tqdm(covg.by_gene(), total=len(covg.genes), desc="Fitting models"):
        model = Model(fpca, fpca_x_values, fpca_basis)
        model.fit(x, var_expl_max, n_pcs_max)
        if model.model is not None:
            models['models'][gene_id] = model
    return models

def fit(batch_covg_dirs: list, batch_id: int, var_expl_max: float, n_pcs_max: int, output_dir: Path, fpca: bool = True, fpca_x_values: str = 'bin', fpca_basis: str = 'discrete'):
    """Fit functional PCA models for all genes in one or more batches"""
    if batch_id is None:
        with open(batch_covg_dirs[0] / 'batches.txt') as f:
            batches = [int(l.strip()) for l in f]
    else:
        batches = [batch_id]
    output_dir.mkdir(exist_ok=True)
    for batch in batches:
        print(f'=== Fitting models for batch {batch} ===', flush=True)
        models = fit_batch(batch_covg_dirs, batch, var_expl_max, n_pcs_max, fpca, fpca_x_values, fpca_basis)
        outfile = output_dir / f'models_{batch}.pickle'
        with open(outfile, 'wb') as f:
            pickle.dump(models, f)
        print(f'Models saved to {outfile}', flush=True)

def transform_batch(batch_covg_dir: Path, batch: int, models_dir: dict) -> Iterator[pd.DataFrame]:
    """Apply functional PCA transformation to all genes in a batch"""
    covg = CoverageData([batch_covg_dir], batch)
    print('  Loading models...', flush=True)
    models_file = models_dir / f'models_{batch}.pickle'
    with open(models_file, 'rb') as f:
        models = pickle.load(f)
    for gene_id, x in tqdm(covg.by_gene(), total=len(covg.genes), desc="  Generating phenotypes per gene"):
        if gene_id in models['models']:
            yield models['models'][gene_id].transform(x)

def transform(batch_covg_dir: Path, models_dir: dict) -> pd.DataFrame:
    """Apply functional PCA transformation to all genes"""
    with open(batch_covg_dir / 'batches.txt') as f:
        batches = [int(l.strip()) for l in f]
    out = []
    for batch in batches:
        print(f'Transforming batch {batch}', flush=True)
        out.extend(transform_batch(batch_covg_dir, batch, models_dir))
    out = pd.concat(out)
    out = out.reset_index()
    return out

def create_parser():
    parser = argparse.ArgumentParser(description='Fit and/or apply models on feature bin coverage data')
    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand', required=True, help='Choose a subcommand')

    # Subparser for 'prepare'
    parser_prepare = subparsers.add_parser('prepare', help='Prepare batched input coverage files. Coverage counts will be assembled and saved in batches, `fit` will be run on each batch separately, and `transform` will run all batches and produce the combined output.')
    parser_prepare.add_argument('-i', '--inputs', type=Path, metavar='FILE', required=True, help='File containing paths to all per-sample coverage files. Each one has one integer per line corresponding to rows of `regions`.')
    parser_prepare.add_argument('-r', '--regions', type=Path, metavar='FILE', required=True, help='BED file containing regions to use for model input data. Must have start, end and region ID in 2nd, 3rd, and 4th columns. Rows must correspond to rows of input coverage files.')
    parser_prepare_phenos = parser_prepare.add_mutually_exclusive_group(required=False)
    parser_prepare_phenos.add_argument('-p', '--pheno-paths', nargs='+', type=Path, metavar='FILE', help='One or more paths to phenotype tables to regress out of the coverage data per gene prior to model input. Files should be in bed format, i.e. input format for tensorqtl. Gene IDs are parsed from the 4th column from the start up to the first non-alphanumeric character.')
    parser_prepare_phenos.add_argument('--pheno-paths-file', type=Path, metavar='FILE', help='File containing list of paths to phenotype tables.')
    parser_prepare.add_argument('-d', '--output-dir', type=Path, metavar='DIR', required=True, help='Directory where per-batch numpy binary files will be written.')
    parser_prepare.add_argument('--batch-size', type=int, default=200, metavar='N', help='Number of genes (at most) per batch. (default: %(default)s)')

    # Subparser for 'fit'
    parser_fit = subparsers.add_parser('fit', help='Fit FPCA or PCA models to coverage data')
    parser_fit_input = parser_fit.add_mutually_exclusive_group(required=True)
    parser_fit_input.add_argument('-d', '--batch-covg-dir', nargs='+', type=Path, metavar='DIR', help='Directory of per-batch numpy binary files. Specify multiple directories to load data from all of them and fit models using the combined dataset. All datasets must have been generated using the same gene bins file.')
    parser_fit_input.add_argument('--dir-file', type=Path, metavar='FILE', help='File containing list of directories of per-batch numpy binary files. Use this instead of -d/--gene-covg-dir in case of many directories.')
    parser_fit.add_argument('-b', '--batch', type=int, metavar='N', help='Batch number to load and fit. Batch numbers start from 0. If omitted, all batches will be loaded and fit in sequence.')
    parser_fit.add_argument('-m', '--models-dir', type=Path, metavar='DIR', required=True, help='Directory in which to save model pickle files.')
    parser_fit.add_argument('-v', '--var-expl-max', type=float, default=0.8, metavar='FLOAT', help='Max variance explained by the PCs kept per gene. Pass 0 or 1 for no variance explained cutoff.  (default: %(default)s)')
    parser_fit.add_argument('-n', '--n-pcs-max', type=int, default=16, metavar='N', help='Max number of PCs to keep per gene. Pass 0 for no cutoff. (default: %(default)s)')
    parser_fit.add_argument('--fpca-x-values', type=str, choices=['bin', 'pos'], default='bin', metavar='STRING', help='Whether to use bin numbers or genomic positions as x-values for functional PCA. (default: %(default)s)')
    parser_fit.add_argument('--fpca-basis', type=str, choices=['discrete', 'spline'], default='discrete', metavar='STRING', help='Basis function to use for functional PCA. `discrete` will run discretized FPCA directly on the data, `spline` will use a 4th-order B-spline basis. (default: %(default)s)')
    parser_fit.add_argument('--regular-pca', action='store_true', help='Use regular PCA instead of functional PCA.')

    # Subparser for 'transform'
    parser_transform = subparsers.add_parser('transform', help='Apply fitted models to coverage data')
    parser_transform.add_argument('-d', '--batch-covg-dir', type=Path, metavar='DIR', required=True, help='Directory of per-batch numpy binary files.')
    parser_transform.add_argument('-m', '--models-dir', type=Path, metavar='DIR', required=True, help='Directory of saved models (*.pickle) to load and use for transformation.')
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
        print('=== Preparing batched coverage files ===', flush=True)
        prepare(args.inputs, args.regions, pheno_paths, args.batch_size, args.output_dir)
        print(f'Coverage files saved in {args.output_dir}', flush=True)
    elif args.subcommand == 'fit':
        if args.batch_covg_dir is not None:
            batch_covg_dirs = args.batch_covg_dir
        else:
            with open(args.dir_file) as f:
                batch_covg_dirs = [Path(l.strip()) for l in f]
        fit(batch_covg_dirs, args.batch, args.var_expl_max, args.n_pcs_max, args.models_dir, not args.regular_pca, args.fpca_x_values, args.fpca_basis)
    elif args.subcommand == 'transform':
        print('=== Generating latent phenotypes ===', flush=True)
        output = transform(args.batch_covg_dir, args.models_dir)
        output.to_csv(args.output, sep='\t', index=False, float_format='%g')
        print(f'Latent phenotypes saved to {args.output}', flush=True)
