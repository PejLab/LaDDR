"""Latent Data-Driven RNA phenotyping"""

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional
import numpy as np
import pandas as pd
import yaml
from .binning import adaptive_binning, fixed_binning, variance_threshold
from .coverage import prepare_coverage
from .init import init_project
from .models import fit, transform, inspect_model
from .setup import setup

@dataclass
class CoverageConfig:
    method: str
    directory: Path
    manifest: Path

@dataclass
class InputConfig:
    gtf: Optional[Path]
    coverage: CoverageConfig
    min_samples_expressed: float
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
    use_existing: bool
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
            max_samples=data_adaptive.get('max_samples', 256),
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
                coverage=coverage,
                min_samples_expressed=data_input.get('min_samples_expressed', 0.5),
                pheno_paths=[Path(p) for p in data_input.get('pheno_paths', [])]
            ),
            binning=BinningConfig(
                use_existing=data_binning.get('use_existing', False),
                method=data_binning.get('method', 'adaptive_diffvar'),
                batch_size=data_binning.get('batch_size', 40),
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
    """Create the CLI parser"""
    parser = argparse.ArgumentParser(description='Fit and/or apply models on feature bin coverage data')
    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand', required=True, help='Choose a subcommand')

    parser_init = subparsers.add_parser('init', help='Initialize a new LaDDR project')
    parser_init.add_argument('project_dir', type=Path, help='Directory to create and initialize project in.')
    parser_init.add_argument('--config-type', type=str, choices=['default', 'example', 'extended'], default='default', help='Type of config file to initialize project with. "default" includes parameters for the recommended binning and model fitting methods. "example" includes a config and coverage manifest file set up to run the example dataset. "extended" includes all possible config parameters.')
    parser_init.add_argument('--template', type=str, choices=['both', 'snakemake', 'shell'], default='both', help='Specify "snakemake" to include a snakefile, "shell" to include a shell script with the basic commands, or "both" to include both.')

    parser_setup = subparsers.add_parser('setup', help='Process annotations and determine gene batches')
    parser_binning = subparsers.add_parser('binning', help='Partition genes into bins for summarizing coverage data')
    parser_prepare = subparsers.add_parser('coverage', help='Prepare RNA-seq coverage for a batch of genes')
    parser_fit = subparsers.add_parser('fit', help='Fit latent phenotype models to coverage data')
    parser_transform = subparsers.add_parser('transform', help='Apply fitted models to coverage data')
    parser_inspect = subparsers.add_parser('inspect', help='Extract model data for visualization and analysis')

    for subparser in [parser_setup, parser_binning, parser_prepare, parser_fit, parser_transform, parser_inspect]:
        subparser.add_argument('-c', '--config', type=Path, metavar='FILE', 
            help='Path to project configuration file. Defaults to config.yaml in project directory.')
        subparser.add_argument('-p', '--project-dir', type=Path, default=Path.cwd(), 
            help='Project directory. Paths in config are relative to this. Defaults to current directory.')

    for subparser in [parser_prepare, parser_transform]:
        subparser.add_argument('-d', '--dataset', type=str, metavar='NAME', help='Name of dataset to process.')

    for subparser in [parser_binning, parser_prepare, parser_fit]:
        subparser.add_argument('-b', '--batch', type=int, metavar='N', help='Batch ID to process. Batch IDs are integers starting from 0. If omitted, all batches will be processed.')

    parser_inspect.add_argument('-g', '--gene-id', type=str, metavar='ID', help='Gene ID of model to inspect.')
    parser_inspect.add_argument('-d', '--dataset', type=str, metavar='NAME', help='Name of dataset whose phenotypes will be used. The output will include mean coverage for the top and bottom tenth of samples for each PC.')

    return parser

def get_sample_table(coverage_config: CoverageConfig, project_dir: Path) -> pd.DataFrame:
    """Get a table of samples from a coverage config"""
    if coverage_config.method == 'directory':
        # Get all subdirectories as datasets, excluding hidden ones
        covg_dir = project_dir / coverage_config.directory
        datasets = [d.name for d in covg_dir.glob('*') 
                   if d.is_dir() and not d.name.startswith('.')]
        # Build table by finding all .bw files in each dataset directory
        rows = []
        for dataset in datasets:
            dataset_dir = covg_dir / dataset
            for bw_file in dataset_dir.glob('*.bw'):
                # Sample ID is filename without .bw extension
                sample = bw_file.stem
                rows.append({
                    'dataset': dataset,
                    'sample': sample,
                    'path': str(bw_file.absolute())
                })
        return pd.DataFrame(rows)
    else:
        manifest_path = project_dir / coverage_config.manifest
        df = pd.read_csv(manifest_path, sep='\t', names=['dataset', 'sample', 'path'])
        duplicates = df['sample'].duplicated(keep=False)
        if duplicates.any():
            duplicate_samples = df[duplicates]
            raise ValueError(
                f"Duplicate sample IDs found in manifest file. Sample IDs must be unique.\n"
                f"Duplicate samples:\n{duplicate_samples.to_string()}"
            )
        # Convert relative paths to absolute paths, preserving existing absolute paths
        if coverage_config.directory:
            # If directory is provided, paths in manifest are relative to that directory
            df['path'] = df['path'].apply(lambda p: str((project_dir / coverage_config.directory / p).absolute()) if not Path(p).is_absolute() else p)
        else:
            # Otherwise paths are relative to project directory
            df['path'] = df['path'].apply(lambda p: str((project_dir / p).absolute()) if not Path(p).is_absolute() else p)
        return df

def cli_setup(config: Config, project_dir: Path, sample_table: pd.DataFrame):
    """Process annotations and determine workflow parameters"""
    assert config.input.gtf is not None, 'gtf is required for setup'
    print('=== Processing annotations, getting batch size, and coverage parameters ===', flush=True)
    genes = setup(
        gtf=project_dir / config.input.gtf,
        bigwig_manifest=sample_table,
        batch_size=config.binning.batch_size,
        outdir=project_dir / 'info'
    )
    if config.binning.method in ['adaptive_covgvar', 'adaptive_diffvar']:
        print('=== Getting binning parameters ===', flush=True)
        # Use coverage from all datasets, subsample if necessary
        with open(project_dir / 'info' / 'median_coverage.txt', 'r') as f:
            median_coverage = float(f.read())
        bigwig_paths = sample_table['path'].to_numpy()
        if len(bigwig_paths) > config.binning.adaptive.max_samples:
            bigwig_paths = np.random.choice(bigwig_paths, config.binning.adaptive.max_samples, replace=False)
        bigwig_paths = [Path(p) for p in bigwig_paths]
        var_per_bin = variance_threshold(
            genes=genes,
            bigwig_paths=bigwig_paths,
            bins_per_gene=config.binning.adaptive.bins_per_gene,
            median_coverage=median_coverage,
            covg_diff=config.binning.method == 'adaptive_diffvar'
        )
        with open(project_dir / 'info' / 'var_per_bin.txt', 'w') as f:
            f.write(f'{var_per_bin:g}')
        print(f"Variance per bin saved to {project_dir / 'info' / 'var_per_bin.txt'}", flush=True)

def cli_binning(args: argparse.Namespace, config: Config, project_dir: Path, sample_table: pd.DataFrame):
    """Partition genes into bins for summarizing coverage data"""
    if config.binning.use_existing:
        raise ValueError('Cannot run binning if use_existing is true')
    print('=== Partitioning genes into bins ===', flush=True)
    if config.binning.method in ['adaptive_covgcorr', 'adaptive_covgvar', 'adaptive_diffvar']:
        # Use coverage from all datasets, subsample if necessary
        bigwig_paths = sample_table['path'].to_numpy()
        if len(bigwig_paths) > config.binning.adaptive.max_samples:
            bigwig_paths = np.random.choice(bigwig_paths, config.binning.adaptive.max_samples, replace=False)
        bigwig_paths = [Path(p) for p in bigwig_paths]
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

    if args.batch is not None:
        batches = [args.batch]
    else:
        with open(project_dir / 'info' / 'n_batches.txt', 'r') as f:
            n_batches = int(f.read())
        batches = list(range(n_batches))

    with open(project_dir / 'info' / 'median_coverage.txt', 'r') as f:
        median_coverage = float(f.read())

    (project_dir / 'covg_norm').mkdir(exist_ok=True)
    
    # Process phenotype paths to replace {dataset} with actual dataset name
    pheno_files = []
    for p in config.input.pheno_paths:
        pheno_path = str(p).format(dataset=dataset)
        pheno_files.append(project_dir / pheno_path)
    
    prepare_coverage(
        bigwig_manifest=sample_table,
        bins_dir=project_dir / 'gene_bins',
        batches=batches,
        pheno_files=pheno_files,
        outdir=project_dir / 'covg_norm' / dataset,
        median_coverage=median_coverage
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
        min_samples_expressed=config.input.min_samples_expressed,
        var_expl_max=config.model.var_explained_max,
        n_pcs_max=config.model.n_pcs_max,
        output_dir=project_dir / 'models',
        fpca=config.model.use_fpca,
        fpca_x_values=config.model.fpca.x_values,
        fpca_basis=config.model.fpca.basis
    )

def cli_transform(args: argparse.Namespace, config: Config, project_dir: Path, sample_table: pd.DataFrame):
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
    output = transform(
        norm_covg_dir=project_dir / 'covg_norm' / dataset,
        models_dir=project_dir / 'models',
        n_batches=n_batches,
        min_samples_expressed=config.input.min_samples_expressed
    )
    outdir = project_dir / 'phenotypes'
    outfile = outdir / f'latent_phenos.{dataset}.tsv.gz'
    outdir.mkdir(exist_ok=True)
    output.to_csv(outfile, sep='\t', index=False, float_format='%g')
    print(f'Latent phenotypes saved to {outfile}', flush=True)

def cli_inspect(args: argparse.Namespace, project_dir: Path, sample_table: pd.DataFrame):
    """Extract model data for visualization and analysis"""
    print(f'=== Extracting model data for {args.gene_id} ===', flush=True)
    if args.dataset is not None:
        dataset = args.dataset
    else:
        datasets = sample_table['dataset'].unique().tolist()
        assert len(datasets) == 1, 'If dataset is omitted, the config must indicate a single dataset'
        dataset = datasets[0]
    print(f'Using coverage and phenotypes from dataset {dataset}', flush=True)
    output = inspect_model(
        gene_id=args.gene_id,
        gene_file=project_dir / 'info' / 'genes.tsv',
        norm_covg_dir=project_dir / 'covg_norm' / dataset,
        models_dir=project_dir / 'models',
        phenotypes=project_dir / 'phenotypes' / f'latent_phenos.{dataset}.tsv.gz',
    )
    outdir = project_dir / 'inspect'
    outfile = outdir / f'inspect.{args.gene_id}.{dataset}.tsv.gz'
    outdir.mkdir(exist_ok=True)
    output.to_csv(outfile, sep='\t', float_format='%g')
    print(f'Model data saved to {outfile}', flush=True)

def cli():
    """LaDDR CLI"""
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
    sample_table = get_sample_table(config.input.coverage, project_dir)

    if args.subcommand == 'setup':
        cli_setup(config, project_dir, sample_table)
    elif args.subcommand == 'binning':
        cli_binning(args, config, project_dir, sample_table)
    elif args.subcommand == 'coverage':
        cli_coverage(args, config, project_dir, sample_table)
    elif args.subcommand == 'fit':
        cli_fit(args, config, project_dir, sample_table)
    elif args.subcommand == 'transform':
        cli_transform(args, config, project_dir, sample_table)
    elif args.subcommand == 'inspect':
        cli_inspect(args, project_dir, sample_table)

if __name__ == '__main__':
    cli()
