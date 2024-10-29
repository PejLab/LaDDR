"""Extract model data for visualization and analysis."""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import pickle
import pyBigWig

class Model:
    pass

def load_coverage(bigwig_paths_file: Path, bins: pd.DataFrame) -> pd.DataFrame:
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
                chrom = str(bins.loc[(gene_id, pos), 'chrom'])
                start = bins.loc[(gene_id, pos), 'chrom_start']
                end = bins.loc[(gene_id, pos), 'chrom_end']
                try:
                    covg[j, i] = f.stats(chrom, start, end, type='mean', exact=True)[0]
                except RuntimeError as e:
                    print(f"RuntimeError for {gene_id}, {chrom}:{start}-{end} in file {bw.stem}: {e}", flush=True)
    samples = [bw.stem for bw in bigwig_paths]
    df = pd.DataFrame(covg, index=bins.index, columns=samples)
    return df

def load_phenos(phenotypes: Path, gene_id: str) -> pd.DataFrame:
    phenos = pd.read_csv(phenotypes, sep='\t', index_col=0)
    phenos = phenos.loc[phenos.index.get_level_values('gene_id') == gene_id]
    phenos.reset_index(level='gene_id', drop=True, inplace=True)
    phenos.set_index('PC', inplace=True)
    phenos = phenos.T.rename_axis('sample_id')
    return phenos

parser = argparse.ArgumentParser(description='Extract loadings from a PCA/fPCA model.')
parser.add_argument('-g', '--gene-id', type=str, metavar='ID', required=True, help='Gene ID to extract.')
parser.add_argument('--bigwig-paths-file', type=Path, metavar='FILE', required=True, help='File containing list of paths to per-sample bigWig files. Basenames of files will be used as sample IDs.')
parser.add_argument('--bins-dir', type=Path, metavar='DIR', required=True, help='Directory of per-batch BED files containing bin regions. Gene batch ID will be retrieved from genes.tsv in this directory.')
parser.add_argument('--norm-covg-dir', type=Path, metavar='DIR', help='Directory of per-batch numpy binary files with normalized coverage. Normalized coverage is not used here, but the accompanying bin info files are used to extract coverage from bigWig files.')
parser.add_argument('--models-dir', type=Path, metavar='DIR', required=True, help='Directory of saved models (*.pickle).')
parser.add_argument('--phenotypes', type=Path, metavar='FILE', required=True, help='Latent phenotype table (TSV). Used to find the top and bottom tenth of samples per PC to get their mean coverage.')
parser.add_argument('-o', '--output', type=Path, metavar='FILE', required=True, help='Output file (TSV).')
args = parser.parse_args()

gene_id = args.gene_id
genes = pd.read_csv(args.bins_dir / 'genes.tsv', sep='\t', index_col=0)
batch = genes.loc[gene_id, 'batch']
bins = pd.read_csv(args.norm_covg_dir / f'batch_{batch}.bins.tsv.gz', sep='\t', index_col=[0, 1])
bins = bins.loc[bins.index.get_level_values('gene_id') == gene_id]
covg = load_coverage(args.bigwig_paths_file, bins)
covg = covg.loc[covg.index.get_level_values('gene_id') == gene_id]
models_file = args.models_dir / f'models_batch_{batch}.pickle'
models = pickle.load(open(models_file, 'rb'))
model = models['models'][gene_id]

if model.fpca:
    # Get the coverage mean of the top and bottom tenth of samples per PC
    phenos = load_phenos(args.phenotypes, gene_id)
    info = covg.iloc[:, :0]
    for pc in phenos.columns:
        top_tenth = phenos[pc].nlargest(int(phenos.shape[0] / 10)).index
        bottom_tenth = phenos[pc].nsmallest(int(phenos.shape[0] / 10)).index
        top_mean = covg.loc[:, top_tenth].mean(axis=1)
        bottom_mean = covg.loc[:, bottom_tenth].mean(axis=1)
        info[f'top_tenth_{pc}'] = top_mean
        info[f'bottom_tenth_{pc}'] = bottom_mean
    if model.fpca_basis == 'discrete':
        loadings = model.model.components_.data_matrix.T
        assert loadings.shape[0] == 1
        loadings = loadings[0]
    else:
        print("Currently unable to extract loadings from fPCA model with basis function.")
        loadings = model.model.components_.coefficients.T
    loadings = pd.DataFrame(loadings, index=info.index, columns=[f'PC{i+1}' for i in range(loadings.shape[1])])
    info = info.join(loadings)
    info.to_csv(args.output, sep='\t')
else:
    # Combine the mean/std from model.features with loadings from the model
    loadings = model.model.components_.T
    loadings = pd.DataFrame(loadings, index=model.features.index, columns=[f'PC{i+1}' for i in range(loadings.shape[1])])
    info = model.features.join(loadings)

    # Get the coverage mean of the top and bottom tenth of samples per PC
    phenos = load_phenos(args.phenotypes, gene_id)
    for pc in phenos.columns:
        top_tenth = phenos[pc].nlargest(int(phenos.shape[0] / 10)).index
        bottom_tenth = phenos[pc].nsmallest(int(phenos.shape[0] / 10)).index
        top_mean = covg.loc[:, top_tenth].mean(axis=1)
        bottom_mean = covg.loc[:, bottom_tenth].mean(axis=1)
        info[f'top_tenth_{pc}'] = top_mean
        info[f'bottom_tenth_{pc}'] = bottom_mean

    info.to_csv(args.output, sep='\t', float_format='%g')
