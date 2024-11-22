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
    method: Optional[str] = 'manifest'
    directory: Optional[Path] = 'covg_bigwig'
    manifest: Optional[Path] = 'coverage_manifest.tsv'

@dataclass
class InputConfig:
    gtf: Optional[Path] = None
    chromosomes: Optional[Path] = None
    coverage: Optional[CoverageConfig] = None
    pheno_paths: Optional[List[Path]] = None

@dataclass 
class AdaptiveBinningConfig:
    max_samples: Optional[int] = 64
    bins_per_gene: Optional[int] = 256
    min_mean_total_covg: Optional[float] = 128
    max_corr: Optional[float] = 0.8

@dataclass
class FixedWidthBinningConfig:
    coding: Optional[int] = 16
    noncoding: Optional[int] = 128

@dataclass
class NbinsPerRegionBinningConfig:
    n_bins: Optional[int] = 24

@dataclass
class BinningConfig:
    method: Optional[str] = 'adaptive3'
    batch_size: Optional[int] = 20
    max_bin_width: Optional[int] = 1024
    adaptive: Optional[AdaptiveBinningConfig] = None
    fixed_width: Optional[FixedWidthBinningConfig] = None
    n_bins_per_region: Optional[NbinsPerRegionBinningConfig] = None

@dataclass
class FPCAConfig:
    x_values: Optional[str] = 'bin'
    basis: Optional[str] = 'discrete'

@dataclass
class ModelConfig:
    var_explained_max: Optional[float] = 0.8
    n_pcs_max: Optional[int] = 16
    use_fpca: Optional[bool] = False
    fpca: Optional[FPCAConfig] = None

@dataclass
class Config:
    input: Optional[InputConfig] = None
    binning: Optional[BinningConfig] = None
    model: Optional[ModelConfig] = None

    @classmethod
    def from_yaml(cls, config_path: Path) -> 'Config':
        with open(config_path) as f:
            data = yaml.safe_load(f)
        # Handle optional top-level configs
        data_input = data.get('input', {})
        data_binning = data.get('binning', {})
        data_model = data.get('model', {})
        return cls(
            input=InputConfig(
                gtf=Path(data_input['gtf']) if 'gtf' in data_input else None,
                chromosomes=Path(data_input['chromosomes']) if 'chromosomes' in data_input else None,
                coverage=CoverageConfig(**data_input.get('coverage', {})),
                pheno_paths=[Path(p) for p in data_input.get('pheno_paths', [])]
            ),
            binning=BinningConfig(
                method=data_binning.get('method'),
                batch_size=data_binning.get('batch_size'),
                max_bin_width=data_binning.get('max_bin_width'),
                adaptive=AdaptiveBinningConfig(**data_binning.get('adaptive', {})),
                fixed_width=FixedWidthBinningConfig(**data_binning.get('fixed_width', {})),
                n_bins_per_region=NbinsPerRegionBinningConfig(**data_binning.get('n_bins_per_region', {}))
            ),
            model=ModelConfig(
                var_explained_max=data_model.get('var_explained_max'),
                n_pcs_max=data_model.get('n_pcs_max'),
                use_fpca=data_model.get('use_fpca'),
                fpca=FPCAConfig(**data_model.get('fpca', {}))
            )
        )

def create_parser():
    parser = argparse.ArgumentParser(description='Fit and/or apply models on feature bin coverage data')
    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand', required=True, help='Choose a subcommand')

    parser_init = subparsers.add_parser('init', help='Initialize a new latent-rna project')
    parser_init.add_argument('project_dir', type=Path, help='Directory to create and initialize project in.')
    parser_init.add_argument('--config-type', type=str, choices=['default', 'example', 'extended'], default='default', help='Type of config file to initialize project with. "default" includes parameters for the recommended binning and model fitting methods. "example" includes a config and coverage manifest file set up to run the example dataset. "extended" includes all possible config parameters.')
    parser_init.add_argument('--template', type=str, choices=['both', 'snakemake', 'shell'], default='both', help='Specify "snakemake" to include a snakefile, "shell" to include a shell script with the basic commands, or "both" to include both.')

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
        init_project(args.project_dir, args.config_type, args.template)
        return
    
    # Look for config file in project directory if not specified
    config_path = args.config if args.config else args.project_dir / 'config.yaml'
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found at {config_path}")
    
    config = Config.from_yaml(config_path)
    project_dir = args.project_dir
    sample_table = get_sample_table(config.input.coverage)

    if args.subcommand == 'binning':
        assert config.input.gtf is not None, 'gtf is required for binning'
        assert config.input.chromosomes is not None, 'chromosomes length file is required for binning'
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
            n_bins_per_region=config.binning.n_bins_per_region.n_bins,
            batch=args.batch
        )

    elif args.subcommand == 'prepare':
        if args.dataset is not None:
            dataset = args.dataset
        else:
            datasets = sample_table['dataset'].unique().tolist()
            assert len(datasets) == 1, 'If dataset is omitted, the config must indicate a single dataset'
            dataset = datasets[0]
        sample_table = sample_table.loc[sample_table['dataset'] == dataset, :]
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
