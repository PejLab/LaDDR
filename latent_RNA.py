"""Run PCA on feature bin coverage and create BED file"""

import argparse
from pathlib import Path
import pickle
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from tqdm import tqdm

class CoverageData:
    def load_sample_files(self, infile: Path, regions: pd.DataFrame, parse_ids: bool):
        self.counts = load_inputs(infile, regions, parse_ids)
        self.samples = self.counts.columns
        self.genes = regions.index.get_level_values('gene_id').unique()

    def from_gene_files(self, gene_covg_dirs: list, regions: pd.DataFrame):
        self.gene_covg_dirs = gene_covg_dirs
        self.counts = None
        self.gene_regions = regions.groupby('gene_id')
        self.samples = []
        for d in self.gene_covg_dirs:
            with open(d / 'samples.txt') as f:
                self.samples.extend([l.strip() for l in f])
        self.genes = self.gene_regions.groups.keys()

    def by_gene(self):
        if self.counts is not None:
            for gene_id, df in self.counts.groupby('gene_id'):
                yield gene_id, df
        else:
            for gene_id in self.genes:
                mats = [np.load(d / f'{gene_id}.npy') for d in self.gene_covg_dirs]
                mat = np.concatenate(mats, axis=1)
                gene_regions = self.gene_regions.get_group(gene_id)
                assert mat.shape[0] == gene_regions.shape[0]
                assert mat.shape[1] == len(self.samples)
                df = pd.DataFrame(mat, index=gene_regions.index, columns=self.samples)
                yield gene_id, df

def load_regions(regionfile: Path) -> pd.DataFrame:
    regions = pd.read_csv(regionfile, sep='\t', header=None, usecols=[1, 2, 3], names=['start', 'end', 'region'], index_col='region')
    regions['length'] = regions['end'] - regions['start']
    regions = regions.drop(['start', 'end'], axis=1)
    regions.index = regions.index.str.split('_', n=1, expand=True)
    regions.index.set_names(['gene_id', 'start'], inplace=True)
    return regions

def load_inputs_np(infiles: list, nrow: int) -> np.ndarray:
    """Load coverage files and concatenate columns into one numpy array"""
    print(f'Loading {len(infiles)} coverage files...', flush=True)
    counts = np.zeros((nrow, len(infiles)), dtype=np.int32)
    for i, fname in tqdm(enumerate(infiles), total=len(infiles), desc="Loading files"):
        d = pd.read_csv(fname, header=None, dtype=np.int32, names=['count'], index_col=False)
        assert d.shape[0] == nrow
        counts[:, i] = d['count']
    return counts

def load_inputs(infile: Path, regions: pd.DataFrame, parse_ids: bool) -> pd.DataFrame:
    """Load coverage files into one dataframe"""
    infiles = pd.read_csv(infile, header=None, names=['path'])
    if parse_ids:
        infiles['sample'] = [Path(p).name.split('.')[0] for p in infiles['path']]
        assert infiles['sample'].is_unique
    else:
        infiles['sample'] = [f'S{i + 1}' for i in range(infiles.shape[0])]
    counts = load_inputs_np(infiles['path'], nrow=regions.shape[0])
    counts = pd.DataFrame(counts, index=regions.index, columns=infiles['sample'])
    return counts

def norm_counts_over_gene(counts: pd.DataFrame) -> pd.DataFrame:
    """For one gene, scale coverage per sample to have the same total per sample
    
    Scale each total to the median of the totals per sample, except for samples
    with 0 total, which are set to 0. This is to prevent NaNs which would cause
    the entire gene to be dropped if it had 0 expression in a single sample.
    """
    sample_totals = counts.sum(axis=0)
    scale_factors = sample_totals.median() / sample_totals
    scale_factors = scale_factors.fillna(0).replace(np.inf, 0)
    counts = counts.mul(scale_factors, axis=1)
    assert not counts.isna().any().any()
    return counts

def normalize_gene(df: pd.DataFrame, regions: pd.DataFrame) -> pd.DataFrame:
    """Normalize coverage for one gene"""
    # Normalize coverage to ignore total expression, i.e. each sample has same total coverage for the gene
    df = norm_counts_over_gene(df)
    # Normalize counts to mean coverage per bp within each region
    df = df.div(regions['length'].loc[df.index], axis=0)
    # Variance stabilization:
    df = np.sqrt(df)
    return df

def prepare(covg: CoverageData, outdir: Path):
    """Split count matrix by gene and save to numpy binary files"""
    outdir.mkdir(exist_ok=True)
    with open(outdir / 'samples.txt', 'w') as f:
        f.write('\n'.join(covg.samples))
    for gene_id, x in tqdm(covg.by_gene(), total=len(covg.genes), desc="Saving per-gene coverage files"):
        counts = x.to_numpy()
        np.save(outdir / f'{gene_id}.npy', counts)

def subset_pcs(pca: PCA, n_pcs_max: int):
    """Subset the PCA model to n_pcs_max PCs to reduce file size"""
    if n_pcs_max != 0 and n_pcs_max < pca.n_components_:
        pca.components_ = pca.components_[:n_pcs_max, :]
        pca.n_components_ = n_pcs_max
        pca.explained_variance_ = pca.explained_variance_[:n_pcs_max]
        pca.explained_variance_ratio_ = pca.explained_variance_ratio_[:n_pcs_max]
        pca.singular_values_ = pca.singular_values_[:n_pcs_max]
        pca.noise_variance_ = None # Could compute if necessary, but omitting otherwise

def gene_fit(df: pd.DataFrame, regions: pd.DataFrame, n_samples_max: int, var_expl: float, n_pcs_max: int) -> PCA:
    """Fit PCA model to normalized coverage for one gene
    
    Saves enough PCs to explain var_expl variance or n_pcs_max, whichever is smaller.
    """
    if n_samples_max != 0 and n_samples_max < df.shape[1]:
        df = df.sample(n=n_samples_max, axis=1)
    df = normalize_gene(df, regions)
    df = df.loc[df.std(axis=1) > 0, :]
    if df.shape[0] == 0:
        return None
    x = df.values.T
    # Center and scale the data:
    xmean = x.mean(axis=0)
    xstd = x.std(axis=0)
    x = (x - xmean) / xstd
    pca = PCA(n_components=var_expl) if var_expl not in {0, 1} else PCA()
    pca.fit(x)
    subset_pcs(pca, n_pcs_max)
    features = pd.DataFrame(index=df.index) # Save feature stats for transforming new data
    features['mean'] = xmean
    features['std'] = xstd
    return {'pca': pca, 'features': features}

def gene_transform(df: pd.DataFrame, regions: pd.DataFrame, pca: PCA, features: pd.DataFrame) -> pd.DataFrame:
    """Apply PCA model to normalized coverage for one gene"""
    df = normalize_gene(df, regions)
    # Filter df to include same features as the model:
    df = df.loc[features.index, :]
    x = df.values.T
    # Center and scale the data using stats from the model-fitting data:
    x = (x - features['mean'].values) / features['std'].values
    mat = pca.transform(x).T
    pc_names = [f'PC{i + 1}' for i in range(mat.shape[0])]
    out = pd.DataFrame(mat, index=pc_names, columns=df.columns)
    out.index.set_names('PC', inplace=True)
    # Add gene_id to index:
    out.reset_index(inplace=True)
    out.insert(0, 'gene_id', df.index.get_level_values('gene_id')[0])
    out.set_index(['gene_id', 'PC'], inplace=True)
    return out

def gene_fit_transform(df: pd.DataFrame, regions: pd.DataFrame, var_expl: float, n_pcs_max: int) -> pd.DataFrame:
    """PCA-transform normalized coverage for one gene.
    
    Saves enough PCs to explain var_expl variance or n_pcs_max, whichever is smaller.
    """
    df = normalize_gene(df, regions)
    df = df.loc[df.std(axis=1) > 0, :]
    if df.shape[0] == 0:
        return None
    x = df.values.T
    # Center and scale the data:
    x = (x - x.mean(axis=0)) / x.std(axis=0)
    pca = PCA(n_components=var_expl) if var_expl not in {0, 1} else PCA()
    mat = pca.fit_transform(x).T
    if n_pcs_max != 0 and n_pcs_max < mat.shape[0]:
        mat = mat[:n_pcs_max, :]
    pc_names = [f'PC{i + 1}' for i in range(mat.shape[0])]
    out = pd.DataFrame(mat, index=pc_names, columns=df.columns)
    out.index.set_names('PC', inplace=True)
    # Add gene_id to index:
    out.reset_index(inplace=True)
    out.insert(0, 'gene_id', df.index.get_level_values('gene_id')[0])
    out.set_index(['gene_id', 'PC'], inplace=True)
    return out

def fit(covg: CoverageData, regions: pd.DataFrame, n_samples_max: int, var_expl_max: float, n_pcs_max: int) -> dict:
    """Fit PCA models for all genes"""
    models = {'var_expl_max': var_expl_max, 'n_pcs_max': n_pcs_max, 'models': {}}
    regions = regions.groupby('gene_id')
    for gene_id, x in tqdm(covg.by_gene(), total=len(covg.genes), desc="Fitting models"):
        reg = regions.get_group(gene_id)
        model = gene_fit(x, reg, n_samples_max, var_expl_max, n_pcs_max)
        if model is not None:
            models['models'][gene_id] = model
    return models

def transform(covg: CoverageData, regions: pd.DataFrame, models: dict) -> pd.DataFrame:
    """Apply PCA transformation to all genes"""
    regions = regions.groupby('gene_id')
    out = []
    for gene_id, x in tqdm(covg.by_gene(), total=len(covg.genes), desc="Applying transformation"):
        if gene_id in models['models']:
            pca = models['models'][gene_id]['pca']
            features = models['models'][gene_id]['features']
            reg = regions.get_group(gene_id)
            out.append(gene_transform(x, reg, pca, features))
    out = pd.concat(out)
    out = out.reset_index()
    return out

def fit_transform(covg: CoverageData, regions: pd.DataFrame, var_expl_max: float, n_pcs_max: int) -> pd.DataFrame:
    """Fit and transform each gene one by one so not all models are stored in memory at once"""
    regions = regions.groupby('gene_id')
    out = []
    for gene_id, x in tqdm(covg.by_gene(), total=len(covg.genes), desc="Fitting and transforming"):
        reg = regions.get_group(gene_id)
        out.append(gene_fit_transform(x, reg, var_expl_max, n_pcs_max))
    out = pd.concat(out)
    out = out.reset_index()
    return out

def create_parser():
    parser = argparse.ArgumentParser(description='Fit and/or apply PCA model on feature bin coverage data')
    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand', required=True, help='Choose a subcommand')

    # Subparser for 'prepare'
    parser_prepare = subparsers.add_parser('prepare', help='Prepare per-gene input coverage files')
    parser_prepare.add_argument('-i', '--inputs', type=Path, metavar='FILE', required=True, help='File containing paths to all input coverage files.')
    parser_prepare.add_argument('-r', '--regions', type=Path, metavar='FILE', required=True, help='BED file containing regions to use for PCA. Must have start, end and region ID in 2nd, 3rd, and 4th columns. Rows must correspond to rows of input coverage files.')
    parser_prepare.add_argument('-d', '--output-dir', type=Path, metavar='DIR', required=True, help='Directory where per-gene numpy binary files will be written.')

    # Subparser for 'fit'
    parser_fit = subparsers.add_parser('fit', help='Fit PCA model')
    parser_fit_input = parser_fit.add_mutually_exclusive_group(required=True)
    parser_fit_input.add_argument('-i', '--inputs', type=Path, metavar='FILE', help='File containing paths to all input coverage files.')
    parser_fit_input.add_argument('-d', '--gene-covg-dir', nargs='+', type=Path, metavar='DIR', help='Directory of per-gene numpy binary files. Specify multiple directories to load them all and treat as a single dataset.')
    parser_fit_input.add_argument('--dir-file', type=Path, metavar='FILE', help='File containing list of directories of per-gene numpy binary files. Use this instead of -d/--gene-covg-dir in case of many directories.')
    parser_fit.add_argument('-r', '--regions', type=Path, metavar='FILE', required=True, help='BED file containing regions to use for PCA. Must have start, end and region ID in 2nd, 3rd, and 4th columns. Rows must correspond to rows of input coverage files.')
    parser_fit.add_argument('--n-samples-max', type=int, default=1024, metavar='N', help='Max number of samples to use for fitting PCA models for efficiency. Used only for preprocessed per-gene inputs. Loaded samples are randomly subsetted if higher than this, done after loading and concatenating data from multiple datasets if applicable, and the sample subsets are chosen independently per gene. Pass 0 for no cutoff. Default 1024.')
    parser_fit.add_argument('-v', '--var-expl-max', type=float, default=0.8, metavar='FLOAT', help='Max variance explained by the PCs kept per gene. Pass 0 or 1 for no variance explained cutoff. Default 0.8.')
    parser_fit.add_argument('-n', '--n-pcs-max', type=int, default=32, metavar='N', help='Max number of PCs to keep per gene. Pass 0 for no cutoff. Default 32.')
    parser_fit.add_argument('-o', '--output', type=Path, metavar='FILE', required=True, help='Output file (*.pickle) to save PCA models')

    # Subparser for 'transform'
    parser_transform = subparsers.add_parser('transform', help='Apply PCA transformation')
    parser_transform_input = parser_transform.add_mutually_exclusive_group(required=True)
    parser_transform_input.add_argument('-i', '--inputs', type=Path, metavar='FILE', help='File containing paths to all input coverage files. Base file name before first "." is sample ID.')
    parser_transform_input.add_argument('-d', '--gene-covg-dir', nargs=1, type=Path, metavar='DIR', help='Directory of per-gene numpy binary files.')
    parser_transform.add_argument('-r', '--regions', type=Path, metavar='FILE', required=True, help='BED file containing regions to use for PCA. Must have start, end and region ID in 2nd, 3rd, and 4th columns. Rows must correspond to rows of input coverage files.')
    parser_transform.add_argument('-m', '--models', type=Path, metavar='FILE', required=True, help='PCA models file (*.pickle) to load and use for transformation')
    parser_transform.add_argument('-o', '--output', type=Path, metavar='FILE', required=True, help='Output file (TSV)')

    # Subparser for 'fit-transform'
    parser_fit_transform = subparsers.add_parser('fit-transform', help='Fit PCA model and apply transformation without saving models')
    parser_fit_transform_input = parser_fit_transform.add_mutually_exclusive_group(required=True)
    parser_fit_transform_input.add_argument('-i', '--inputs', type=Path, metavar='FILE', help='File containing paths to all input coverage files. Base file name before first "." is sample ID.')
    parser_fit_transform_input.add_argument('-d', '--gene-covg-dir', nargs=1, type=Path, metavar='DIR', help='Directory of per-gene numpy binary files.')
    parser_fit_transform.add_argument('-r', '--regions', type=Path, metavar='FILE', required=True, help='BED file containing regions to use for PCA. Must have start, end and region ID in 2nd, 3rd, and 4th columns. Rows must correspond to rows of input coverage files.')
    parser_fit_transform.add_argument('-v', '--var-expl-max', type=float, default=0.8, metavar='FLOAT', help='Max variance explained by the PCs kept per gene. Pass 0 or 1 for no variance explained cutoff. Default 0.8.')
    parser_fit_transform.add_argument('-n', '--n-pcs-max', type=int, default=32, metavar='N', help='Max number of PCs to keep per gene. Pass 0 for no cutoff. Default 32.')
    parser_fit_transform.add_argument('-o', '--output', type=Path, metavar='FILE', required=True, help='Output file (TSV)')

    return parser

if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()
    regions = load_regions(args.regions)
    # inputs = load_inputs(args.inputs, regions, parse_ids=('transform' in args.subcommand))
    # inputs = CoverageData(regions, args.inputs, args.gene_covg_dir, parse_ids=('transform' in args.subcommand))
    covg = CoverageData()
    if args.subcommand == 'prepare':
        covg.load_sample_files(args.inputs, regions, parse_ids=True)
        prepare(covg, args.output_dir)
    else:
        if args.inputs is not None:
            covg.load_sample_files(args.inputs, regions, parse_ids=('transform' in args.subcommand))
        elif args.gene_covg_dir is not None:
            covg.from_gene_files(args.gene_covg_dir, regions)
        else:
            with open(args.dir_file) as f:
                gene_covg_dirs = [Path(l.strip()) for l in f]
            covg.from_gene_files(gene_covg_dirs, regions)
    if args.subcommand == 'fit':
        print('Fitting PCA models...', flush=True)
        models = fit(covg, regions, args.n_samples_max, args.var_expl_max, args.n_pcs_max)
        with open(args.output, 'wb') as f:
            pickle.dump(models, f)
        print(f'PCA models saved to {args.output}', flush=True)
    elif args.subcommand == 'transform':
        print('Loading PCA models...', flush=True)
        with open(args.models, 'rb') as f:
            models = pickle.load(f)
        print('Applying PCA transformation...', flush=True)
        output = transform(covg, regions, models)
        output.to_csv(args.output, sep='\t', index=False, float_format='%g')
        print(f'PCA-transformed features saved to {args.output}', flush=True)
    elif args.subcommand == 'fit-transform':
        print('Fitting PCA models and applying transformation...', flush=True)
        output = fit_transform(covg, regions, args.var_expl_max, args.n_pcs_max)
        output.to_csv(args.output, sep='\t', index=False, float_format='%g')
        print(f'PCA-transformed features saved to {args.output}', flush=True)
