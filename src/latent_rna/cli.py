"""Extract latent transcriptomic phenotypes from RNA-seq coverage data"""

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import List
import numpy as np
import pandas as pd
import yaml
from .binning import adaptive_binning, fixed_binning, variance_threshold
from .coverage import prepare_coverage
from .init import init_project
from .models import fit, transform
from .setup import setup, compute_sample_scaling_factors

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
    pheno_paths: List[Path]

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
class FixedCountBinningConfig:
    n_bins_per_region: int

@dataclass
class BinningConfig:
    method: str
    batch_size: int
    max_bin_width: int
    adaptive: AdaptiveBinningConfig
    fixed_width: FixedWidthBinningConfig
    fixed_count: FixedCountBinningConfig

@dataclass
class FPCAConfig:
    x_values: str
    basis: str

@dataclass
class ModelConfig:
    use_existing: bool
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
        # Handle optional top-level configs
        data_input = data.get('input', {})
        data_binning = data.get('binning', {})
        data_model = data.get('model', {})
        # Handle optional nested configs
        data_coverage = data_input.get('coverage', {})
        coverage = CoverageConfig(
            method=data_coverage.get('method', 'manifest'),
            directory=Path(data_coverage.get('directory', 'covg_bigwig')),
            manifest=Path(data_coverage.get('manifest', 'coverage_manifest.tsv'))
        )
        data_adaptive = data_binning.get('adaptive', {})
        adaptive = AdaptiveBinningConfig(
            max_samples=data_adaptive.get('max_samples', 64),
            bins_per_gene=data_adaptive.get('bins_per_gene', 256),
            min_mean_total_covg=data_adaptive.get('min_mean_total_covg', 128),
            max_corr=data_adaptive.get('max_corr', 0.8)
        )
        data_fixed_width = data_binning.get('fixed_width', {})
        fixed_width = FixedWidthBinningConfig(
            coding=data_fixed_width.get('coding', 16),
            noncoding=data_fixed_width.get('noncoding', 128)
        )
        data_fixed_count = data_binning.get('fixed_count', {})
        fixed_count = FixedCountBinningConfig(
            n_bins_per_region=data_fixed_count.get('n_bins_per_region', 24)
        )
        data_fpca = data_model.get('fpca', {})
        fpca = FPCAConfig(
            x_values=data_fpca.get('x_values', 'bin'),
            basis=data_fpca.get('basis', 'discrete')
        )
        return cls(
            input=InputConfig(
                gtf=Path(data_input['gtf']) if 'gtf' in data_input else None,
                chromosomes=Path(data_input['chromosomes']) if 'chromosomes' in data_input else None,
                coverage=coverage,
                pheno_paths=[Path(p) for p in data_input.get('pheno_paths', [])]
            ),
            binning=BinningConfig(
                method=data_binning.get('method', 'adaptive_diffvar'),
                batch_size=data_binning.get('batch_size', 20),
                max_bin_width=data_binning.get('max_bin_width', 1024),
                adaptive=adaptive,
                fixed_width=fixed_width,
                fixed_count=fixed_count
            ),
            model=ModelConfig(
                use_existing=data_model.get('use_existing', False),
                var_explained_max=data_model.get('var_explained_max', 0.8),
                n_pcs_max=data_model.get('n_pcs_max', 16),
                use_fpca=data_model.get('use_fpca', False),
                fpca=fpca
            )
        )

def create_parser():
    parser = argparse.ArgumentParser(description='Fit and/or apply models on feature bin coverage data')
    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand', required=True, help='Choose a subcommand')

    parser_init = subparsers.add_parser('init', help='Initialize a new latent-rna project')
    parser_init.add_argument('project_dir', type=Path, help='Directory to create and initialize project in.')
    parser_init.add_argument('--config-type', type=str, choices=['default', 'example', 'extended'], default='default', help='Type of config file to initialize project with. "default" includes parameters for the recommended binning and model fitting methods. "example" includes a config and coverage manifest file set up to run the example dataset. "extended" includes all possible config parameters.')
    parser_init.add_argument('--template', type=str, choices=['both', 'snakemake', 'shell'], default='both', help='Specify "snakemake" to include a snakefile, "shell" to include a shell script with the basic commands, or "both" to include both.')

    parser_setup = subparsers.add_parser('setup', help='Process annotations and determine gene batches')
    parser_binning = subparsers.add_parser('binning', help='Partition genes into bins for summarizing coverage data')
    parser_prepare = subparsers.add_parser('coverage', help='Prepare RNA-seq coverage for a batch of genes')
    parser_fit = subparsers.add_parser('fit', help='Fit latent phenotype models to coverage data')
    parser_transform = subparsers.add_parser('transform', help='Apply fitted models to coverage data')

    for subparser in [parser_setup, parser_binning, parser_prepare, parser_fit, parser_transform]:
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

def cli_setup(config: Config, project_dir: Path, sample_table: pd.DataFrame):
    """Process annotations and determine gene batches"""
    assert config.input.gtf is not None, 'gtf is required for setup'
    assert config.input.chromosomes is not None, 'chromosomes length file is required for setup'
    genes = setup(
        gtf=project_dir / config.input.gtf,
        chrom_file=project_dir / config.input.chromosomes,
        batch_size=config.binning.batch_size,
        outdir=project_dir / 'info'
    )
    if config.binning.method in ['adaptive_covgvar', 'adaptive_diffvar']:
        # Use coverage from all datasets, subsample if necessary
        bigwig_paths = [Path(p) for p in sample_table['path'].tolist()]
        if len(bigwig_paths) > config.binning.adaptive.max_samples:
            bigwig_paths = np.random.choice(bigwig_paths, config.binning.adaptive.max_samples, replace=False)
        var_per_bin = variance_threshold(
            genes=genes,
            bigwig_paths=bigwig_paths,
            bins_per_gene=config.binning.adaptive.bins_per_gene,
            covg_diff=config.binning.method == 'adaptive_diffvar'
        )
        with open(project_dir / 'info' / 'var_per_bin.txt', 'w') as f:
            f.write(f'{var_per_bin:g}')
    compute_sample_scaling_factors(
        bigwig_manifest=sample_table,
        genes=genes,
        outdir=project_dir / 'info',
        use_existing=config.model.use_existing
    )

def cli_binning(args: argparse.Namespace, config: Config, project_dir: Path, sample_table: pd.DataFrame):
    """Partition genes into bins for summarizing coverage data"""
    if config.binning.method in ['adaptive_covgcorr', 'adaptive_covgvar', 'adaptive_diffvar']:
        # Use coverage from all datasets, subsample if necessary
        bigwig_paths = [Path(p) for p in sample_table['path'].tolist()]
        if len(bigwig_paths) > config.binning.adaptive.max_samples:
            bigwig_paths = np.random.choice(bigwig_paths, config.binning.adaptive.max_samples, replace=False)
        if config.binning.method in {'adaptive_covgvar', 'adaptive_diffvar'}:
            with open(project_dir / 'info' / 'var_per_bin.txt', 'r') as f:
                var_threshold = float(f.read())
        else:
            var_threshold = None
        adaptive_binning(
            gene_file=project_dir / 'info' / 'genes.tsv',
            exon_file=project_dir / 'info' / 'exons.tsv.gz',
            outdir=project_dir / 'gene_bins',
            binning_method=config.binning.method,
            bigwig_paths=bigwig_paths,
            var_threshold=var_threshold,
            min_mean_total_covg=config.binning.adaptive.min_mean_total_covg,
            max_corr=config.binning.adaptive.max_corr,
            max_bin_width=config.binning.max_bin_width,
            batch=args.batch
        )
    elif config.binning.method in ['fixed_width', 'fixed_count']:
        fixed_binning(
            gene_file=project_dir / 'info' / 'genes.tsv',
            exon_file=project_dir / 'info' / 'exons.tsv.gz',
            outdir=project_dir / 'gene_bins',
            binning_method=config.binning.method,
            max_bin_width=config.binning.max_bin_width,
            bin_width_coding=config.binning.fixed_width.coding,
            bin_width_noncoding=config.binning.fixed_width.noncoding,
            n_bins_per_region=config.binning.fixed_count.n_bins_per_region,
            batch=args.batch
        )
    else:
        raise ValueError(f'Invalid binning method: {config.binning.method}')

def cli_coverage(args: argparse.Namespace, config: Config, project_dir: Path, sample_table: pd.DataFrame):
    """Prepare RNA-seq coverage for a batch of genes"""
    if args.dataset is not None:
        dataset = args.dataset
    else:
        datasets = sample_table['dataset'].unique().tolist()
        assert len(datasets) == 1, 'If dataset is omitted, the config must indicate a single dataset'
        dataset = datasets[0]
    sample_table = sample_table.loc[sample_table['dataset'] == dataset, :]

    scaling_factors = pd.read_csv(project_dir / 'info' / 'scaling_factors.tsv', sep='\t', index_col='sample')['scaling_factor']
    missing_samples = set(sample_table['sample']) - set(scaling_factors.index)
    if missing_samples:
        raise ValueError(
            f"The following samples are missing scaling factors: {missing_samples}\n"
            "This may occur if setup was run with a different set of samples.\n"
            "If applying existing models to new samples, set model.use_existing to true in the config."
        )
    if args.batch is not None:
        batches = [args.batch]
    else:
        with open(project_dir / 'info' / 'n_batches.txt', 'r') as f:
            n_batches = int(f.read())
        batches = list(range(n_batches))
    (project_dir / 'covg_norm').mkdir(exist_ok=True)
    prepare_coverage(
        bigwig_manifest=sample_table,
        scaling_factors=scaling_factors,
        bins_dir=project_dir / 'gene_bins',
        batches=batches,
        pheno_files=[project_dir / p for p in config.input.pheno_paths],
        outdir=project_dir / 'covg_norm' / dataset
    )

def cli_fit(args: argparse.Namespace, config: Config, project_dir: Path, sample_table: pd.DataFrame):
    """Fit latent phenotype models to normalized coverage data"""
    if config.model.use_existing:
        raise ValueError('Cannot fit models if use_existing is true')
    datasets = sample_table['dataset'].unique().tolist()
    if args.batch is not None:
        batches = [args.batch]
    else:
        with open(project_dir / 'info' / 'n_batches.txt', 'r') as f:
            n_batches = int(f.read())
        batches = list(range(n_batches))
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

def cli_transform(args: argparse.Namespace, project_dir: Path, sample_table: pd.DataFrame):
    """Apply fitted models to normalized coverage data"""
    if args.dataset is not None:
        dataset = args.dataset
    else:
        datasets = sample_table['dataset'].unique().tolist()
        assert len(datasets) == 1, 'If dataset is omitted, the config must indicate a single dataset'
        dataset = datasets[0]
    with open(project_dir / 'info' / 'n_batches.txt', 'r') as f:
        n_batches = int(f.read())
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

def cli():
    """Latent RNA CLI"""
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

    if args.subcommand == 'setup':
        cli_setup(config, project_dir, sample_table)

    elif args.subcommand == 'binning':
        cli_binning(args, config, project_dir, sample_table)

    elif args.subcommand == 'coverage':
        cli_coverage(args, config, project_dir, sample_table)
    elif args.subcommand == 'fit':
        cli_fit(args, config, project_dir, sample_table)

    elif args.subcommand == 'transform':
        cli_transform(args, project_dir, sample_table)

if __name__ == '__main__':
    cli()
