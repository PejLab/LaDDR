"""Extract latent transcriptomic phenotypes from RNA-seq coverage data"""

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List
import numpy as np
import pandas as pd
import yaml
from .binning import binning
from .coverage import prepare
from .init import init_project
from .models import fit, transform

@dataclass
class CoverageConfig:
    method: str
    directory: Path
    manifest: Path

@dataclass
class InputConfig:
    gtf: Path
    chromosomes: Path
    coverage: CoverageConfig
    pheno_paths: Optional[List[Path]] = None

@dataclass 
class AdaptiveBinningConfig:
    max_samples: int
    bins_per_gene: int
    min_mean_total_covg: float
    max_corr: float

@dataclass
class FixedWidthBinningConfig:
    coding: int
    noncoding: int

@dataclass
class BinningConfig:
    method: str
    batch_size: int
    max_bin_width: int
    n_bins: int
    adaptive: AdaptiveBinningConfig
    fixed_width: FixedWidthBinningConfig

@dataclass
class FPCAConfig:
    x_values: str
    basis: str

@dataclass
class ModelConfig:
    var_explained_max: float
    n_pcs_max: int
    use_fpca: bool
    fpca: FPCAConfig

@dataclass
class Config:
    input: InputConfig
    binning: BinningConfig
    model: ModelConfig

    @classmethod
    def from_yaml(cls, config_path: Path) -> 'Config':
        with open(config_path) as f:
            data = yaml.safe_load(f)
        return cls(
            input=InputConfig(
                gtf=Path(data['input']['gtf']),
                chromosomes=Path(data['input']['chromosomes']),
                coverage=CoverageConfig(**data['input']['coverage']),
                pheno_paths=[Path(p) for p in data['input']['pheno_paths']] 
                    if data['input']['pheno_paths'] else []
            ),
            binning=BinningConfig(
                method=data['binning']['method'],
                batch_size=data['binning']['batch_size'],
                max_bin_width=data['binning']['max_bin_width'],
                n_bins=data['binning']['n_bins'],
                adaptive=AdaptiveBinningConfig(**data['binning']['adaptive']),
                fixed_width=FixedWidthBinningConfig(**data['binning']['fixed_width'])
            ),
            model=ModelConfig(
                var_explained_max=data['model']['var_explained_max'],
                n_pcs_max=data['model']['n_pcs_max'],
                use_fpca=data['model']['use_fpca'],
                fpca=FPCAConfig(**data['model']['fpca'])
            )
        )

def create_parser():
    parser = argparse.ArgumentParser(description='Fit and/or apply models on feature bin coverage data')
    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand', required=True, help='Choose a subcommand')

    parser_init = subparsers.add_parser('init', help='Initialize a new latent-rna project')
    parser_init.add_argument('project_dir', type=Path, help='Directory to create and initialize project in.')

    parser_binning = subparsers.add_parser('binning', help='Partition genes into bins for summarizing coverage data')
    parser_prepare = subparsers.add_parser('prepare', help='Prepare RNA-seq coverage for a batch of genes')
    parser_fit = subparsers.add_parser('fit', help='Fit FPCA or PCA models to coverage data')
    parser_transform = subparsers.add_parser('transform', help='Apply fitted models to coverage data')

    for subparser in [parser_binning, parser_prepare, parser_fit, parser_transform]:
        subparser.add_argument('-c', '--config', type=Path, metavar='FILE', 
            help='Path to project configuration file. Defaults to config.yaml in project directory.')
        subparser.add_argument('-p', '--project-dir', type=Path, default=Path.cwd(), 
            help='Project directory. Paths in config are relative to this. Defaults to current directory.')

    for subparser in [parser_prepare, parser_transform]:
        subparser.add_argument('-d', '--dataset', type=str, metavar='NAME', help='Name of dataset to process.')

    for subparser in [parser_binning, parser_prepare, parser_fit]:
        subparser.add_argument('-b', '--batch', type=int, metavar='N', help='Batch ID to process. Batch IDs are integers starting from 0. If omitted, all batches will be processed.')

    return parser

def get_sample_table(coverage_config: CoverageConfig) -> pd.DataFrame:
    if coverage_config.method == 'directory':
        # Get all subdirectories as datasets, excluding hidden ones
        datasets = [d.name for d in coverage_config.directory.glob('*') 
                   if d.is_dir() and not d.name.startswith('.')]
        # Build table by finding all .bw files in each dataset directory
        rows = []
        for dataset in datasets:
            dataset_dir = coverage_config.directory / dataset
            for bw_file in dataset_dir.glob('*.bw'):
                # Sample ID is filename without .bw extension
                sample = bw_file.stem
                rows.append({
                    'dataset': dataset,
                    'sample': sample,
                    'path': str(bw_file)
                })
        return pd.DataFrame(rows)
    else:
        return pd.read_csv(coverage_config.manifest, sep='\t', names=['dataset', 'sample', 'path'])

def cli():
    parser = create_parser()
    args = parser.parse_args()
    if args.subcommand == 'init':
        init_project(args.project_dir)
        return
    
    # Look for config file in project directory if not specified
    config_path = args.config if args.config else args.project_dir / 'config.yaml'
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found at {config_path}")
    
    config = Config.from_yaml(config_path)
    project_dir = args.project_dir
    sample_table = get_sample_table(config.input.coverage)

    if args.subcommand == 'binning':
        # Use coverage from all datasets, subsample if necessary
        bigwig_paths = sample_table['path'].tolist()
        if config.binning.adaptive.max_samples and len(bigwig_paths) > config.binning.adaptive.max_samples:
            bigwig_paths = np.random.choice(bigwig_paths, config.binning.adaptive.max_samples, replace=False)
        binning(
            gtf=project_dir / config.input.gtf,
            chromosomes=project_dir / config.input.chromosomes,
            outdir=project_dir / 'gene_bins',
            batch_size=config.binning.batch_size,
            binning_method=config.binning.method,
            bigwig_paths=bigwig_paths,
            bins_per_gene=config.binning.adaptive.bins_per_gene,
            min_mean_total_covg=config.binning.adaptive.min_mean_total_covg,
            max_corr=config.binning.adaptive.max_corr,
            max_bin_width=config.binning.max_bin_width,
            bin_width_coding=config.binning.fixed_width.coding,
            bin_width_noncoding=config.binning.fixed_width.noncoding,
            n_bins=config.binning.n_bins,
            batch=args.batch
        )

    elif args.subcommand == 'prepare':
        if args.dataset is not None:
            dataset = args.dataset
        else:
            datasets = sample_table['dataset'].unique().tolist()
            assert len(datasets) == 1, 'If dataset is omitted, the config must indicate a single dataset'
            dataset = datasets[0]
        n_batches = len(list(project_dir.glob('gene_bins/batch_*.bed.gz')))
        batches = [args.batch] if args.batch is not None else list(range(n_batches))
        (project_dir / 'covg_norm').mkdir(exist_ok=True)
        prepare(
            bigwig_manifest=sample_table,
            bins_dir=project_dir / 'gene_bins',
            batches=batches,
            pheno_files=[project_dir / p for p in config.input.pheno_paths],
            outdir=project_dir / 'covg_norm' / dataset
        )

    elif args.subcommand == 'fit':
        datasets = sample_table['dataset'].unique().tolist()
        n_batches = len(list(project_dir.glob('gene_bins/batch_*.bed.gz')))
        batches = [args.batch] if args.batch is not None else list(range(n_batches))
        fit(
            norm_covg_dirs=[project_dir / 'covg_norm' / dataset for dataset in datasets],
            batches=batches,
            var_expl_max=config.model.var_explained_max,
            n_pcs_max=config.model.n_pcs_max,
            output_dir=project_dir / 'models',
            fpca=config.model.use_fpca,
            fpca_x_values=config.model.fpca.x_values,
            fpca_basis=config.model.fpca.basis
        )

    elif args.subcommand == 'transform':
        if args.dataset is not None:
            dataset = args.dataset
        else:
            datasets = sample_table['dataset'].unique().tolist()
            assert len(datasets) == 1, 'If dataset is omitted, the config must indicate a single dataset'
            dataset = datasets[0]
        n_batches = len(list(project_dir.glob('models/models_batch_*.pickle')))
        print('=== Generating latent phenotypes ===', flush=True)
        outdir = project_dir / 'phenotypes'
        outfile = outdir / f'latent_phenos.{dataset}.tsv.gz'
        output = transform(
            norm_covg_dir=project_dir / 'covg_norm' / dataset,
            models_dir=project_dir / 'models',
            n_batches=n_batches
        )
        outdir.mkdir(exist_ok=True)
        output.to_csv(outfile, sep='\t', index=False, float_format='%g')
        print(f'Latent phenotypes saved to {outfile}', flush=True)

if __name__ == '__main__':
    cli()
