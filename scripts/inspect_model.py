"""Extract model data for visualization and analysis."""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import pickle

class Model:
    pass

def load_coverage(batch_covg_dir: Path, batch_id: str) -> pd.DataFrame:
    mat = np.load(batch_covg_dir / f'covg_{batch_id}.npy')
    region_file = batch_covg_dir / f'covg_{batch_id}.regions.tsv.gz'
    regions = pd.read_csv(region_file, sep='\t', index_col=[0, 1])
    with open(batch_covg_dir / 'samples.txt', 'r') as f:
        samples = f.read().splitlines()
    assert mat.shape[0] == regions.shape[0]
    assert mat.shape[1] == len(samples)
    df = pd.DataFrame(mat, index=regions.index, columns=samples)
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
parser.add_argument('-b', '--batch', type=str, metavar='ID', required=True, help='Batch ID containing the gene.')
parser.add_argument('-c', '--batch-covg-dir', type=Path, metavar='DIR', required=True, help='Directory of per-batch numpy binary files.')
parser.add_argument('-m', '--models-dir', type=Path, metavar='DIR', required=True, help='Directory of saved models (*.pickle).')
parser.add_argument('-p', '--phenotypes', type=Path, metavar='FILE', required=True, help='Latent phenotype table (TSV). Used to find the top and bottom tenth of samples per PC to get their mean coverage.')
parser.add_argument('-o', '--output', type=Path, metavar='FILE', required=True, help='Output file (TSV).')
args = parser.parse_args()

gene_id = args.gene_id
covg = load_coverage(args.batch_covg_dir, args.batch)
covg = covg.loc[covg.index.get_level_values('gene_id') == gene_id]
models_file = args.models_dir / f'models_{args.batch}.pickle'
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

    info.to_csv(args.output, sep='\t')
