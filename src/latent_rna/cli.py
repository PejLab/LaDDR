"""Extract latent transcriptomic phenotypes from RNA-seq coverage data"""

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List
import yaml
from .binning import binning
from .coverage import prepare
from .init import init_project
from .models import fit, transform

@dataclass
class InputConfig:
    gtf: Path
    chromosomes: Path
    bigwig_paths: Path
    pheno_paths: Optional[List[Path]] = None

@dataclass 
class AdaptiveBinningConfig:
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
class ModelConfig:
    var_explained_max: float
    n_pcs_max: int
    use_fpca: bool
    fpca_x_values: str = 'bin'
    fpca_basis: str = 'discrete'

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
                bigwig_paths=Path(data['input']['bigwig_paths']),
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
            model=ModelConfig(**data['model'])
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
        subparser.add_argument('-c', '--config', type=Path, required=True, metavar='FILE', help='Path to project configuration file.')
        subparser.add_argument('-p', '--project-dir', type=Path, default=Path.cwd(), help='Project directory. Paths in config are relative to this. Defaults to current directory.')

    for subparser in [parser_binning, parser_prepare, parser_transform]:
        subparser.add_argument('-d', '--dataset', type=Path, metavar='FILE', help='Name of dataset to process.')

    for subparser in [parser_binning, parser_prepare, parser_fit]:
        subparser.add_argument('-b', '--batch', type=int, metavar='N', help='Batch ID to process. Batch IDs are integers starting from 0. If omitted, all batches will be processed.')

    return parser

def cli():
    parser = create_parser()
    args = parser.parse_args()
    if args.subcommand == 'init':
        init_project(args.project_dir)
        return
    config = Config.from_yaml(args.config)
    project_dir = args.project_dir
    if args.subcommand == 'binning':
        if args.dataset is not None:
            dataset = args.dataset
        else:
            # Assert project_dir / covg_bigwig contains only one subdirectory and use that as dataset
            assert len(list(project_dir.glob('covg_bigwig/*'))) == 1, 'If dataset is omitted, there should be exactly one subdirectory in project_dir/covg_bigwig'
            dataset = next(project_dir.glob('covg_bigwig/*')).name
        binning(
            gtf=project_dir / config.input.gtf,
            chromosomes=project_dir / config.input.chromosomes,
            outdir=project_dir / 'gene_bins' / dataset,
            batch_size=config.binning.batch_size,
            binning_method=config.binning.method,
            bigwig_paths_file=project_dir / config.input.bigwig_paths,
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
            # Assert project_dir / covg_bigwig contains only one subdirectory and use that as dataset
            assert len(list(project_dir.glob('covg_bigwig/*'))) == 1, 'If dataset is omitted, there should be exactly one subdirectory in project_dir/covg_bigwig'
            dataset = next(project_dir.glob('covg_bigwig/*')).name
        batches = [args.batch] if args.batch is not None else list(range(args.n_batches))
        prepare(
            bigwig_paths_file=project_dir / config.input.bigwig_paths,
            bins_dir=project_dir / 'gene_bins' / dataset,
            batches=batches,
            pheno_files=[project_dir / p for p in config.input.pheno_paths],
            outdir=project_dir / 'covg_norm' / dataset
        )
    elif args.subcommand == 'fit':
        batches = [args.batch] if args.batch is not None else list(range(args.n_batches))
        datasets = list(project_dir.glob('gene_bins/*'))
        fit(
            norm_covg_dirs=[project_dir / 'covg_norm' / p for p in datasets],
            batches=batches,
            var_expl_max=config.model.var_explained_max,
            n_pcs_max=config.model.n_pcs_max,
            output_dir=project_dir / 'models',
            fpca=config.model.use_fpca,
            fpca_x_values=config.model.fpca_x_values,
            fpca_basis=config.model.fpca_basis
        )
    elif args.subcommand == 'transform':
        if args.dataset is not None:
            dataset = args.dataset
        else:
            # Assert project_dir / covg_bigwig contains only one subdirectory and use that as dataset
            assert len(list(project_dir.glob('covg_norm/*'))) == 1, 'If dataset is omitted, there should be exactly one subdirectory in project_dir/covg_norm'
            dataset = next(project_dir.glob('covg_norm/*')).name
        n_batches = len(list(project_dir.glob('models/models_batch_*.pickle')))
        print('=== Generating latent phenotypes ===', flush=True)
        outfile = project_dir / 'phenotypes' / f'latent_phenotypes_{dataset}.tsv.gz'
        output = transform(
            norm_covg_dir=project_dir / 'covg_norm' / dataset,
            models_dir=project_dir / 'models',
            n_batches=n_batches
        )
        output.to_csv(outfile, sep='\t', index=False, float_format='%g')
        print(f'Latent phenotypes saved to {outfile}', flush=True)

if __name__ == '__main__':
    cli()
