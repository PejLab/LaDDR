"""Extract latent transcriptomic phenotypes from RNA-seq coverage data"""

import argparse
from pathlib import Path
from .binning import binning
from .coverage import prepare
from .models import fit, transform

def create_parser():
    parser = argparse.ArgumentParser(description='Fit and/or apply models on feature bin coverage data')
    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand', required=True, help='Choose a subcommand')

    # Subparser for 'binning'
    parser_binning = subparsers.add_parser('binning', help='Partition genes into bins for summarizing coverage data')
    parser_binning.add_argument('-g', '--gtf', type=Path, required=True, metavar='FILE', help='Transcript annotation in GTF format. Must be the collapsed annotation produced by `collapse_annotation.py`.')
    parser_binning.add_argument('-c', '--chromosomes', type=Path, required=True, metavar='FILE', help='Chromosome lengths file, e.g. chrNameLength.txt from STAR index, to sort chromosomes.')
    parser_binning.add_argument('--outdir', type=Path, required=True, metavar='PATH', help='Directory in which to save per-batch BED files.')
    parser_binning.add_argument('--binning-method', choices=['adaptive1', 'adaptive2', 'adaptive3', 'bin-width', 'n-bins'], default='n-bins', help='Whether to determine bins adaptively from coverage data, split all coding/noncoding regions into fixed-width bins, or split each region into a fixed number of bins. (default: %(default)s)')
    parser_binning.add_argument('--bigwig-paths-file', type=Path, metavar='FILE', help='For "adaptive" methods, file containing list of paths to per-sample bigWig files to use for coverage data.')
    parser_binning.add_argument('--min-mean-total-covg', type=float, default=128, metavar='FLOAT', help='For method "adaptive1", minimum allowed mean total coverage per sample for a bin. (default: %(default)s)')
    parser_binning.add_argument('--max-corr', type=float, default=0.8, metavar='FLOAT', help='For method "adaptive1", maximum allowed correlation between normalized coverage of adjacent bins. (default: %(default)s)')
    parser_binning.add_argument('--bins-per-gene', type=int, default=128, metavar='N', help='For method "adaptive2" or "adaptive3", approximate number of bins to create per gene on average. (default: %(default)s)')
    parser_binning.add_argument('--bin-width-coding', type=int, default=16, metavar='N', help='For method "bin-width", width of bins for coding (exonic) regions of the genes. (default: %(default)s)')
    parser_binning.add_argument('--bin-width-noncoding', type=int, default=128, metavar='N', help='For method "bin-width", width of bins for noncoding (non-exonic) regions of the genes. (default: %(default)s)')
    parser_binning.add_argument('--n-bins', type=int, default=24, metavar='N', help='For method "n-bins", number of bins to split each feature (exon, intron, etc.) into. (default: %(default)s)')
    parser_binning.add_argument('--max-bin-width', type=int, default=1024, metavar='N', help='After a binning method is run, any bins larger than this will be split up. (default: %(default)s)')
    parser_binning.add_argument('--batch-size', type=int, default=200, metavar='N', help='Number of genes (at most) per batch. (default: %(default)s)')
    parser_binning.add_argument('--batch', type=int, metavar='N', help='Batch ID to process. Batch IDs are integers starting from 0. If omitted, all batches will be processed.')

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

def cli():
    parser = create_parser()
    args = parser.parse_args()
    if args.subcommand == 'binning':
        binning(args.gtf, args.chromosomes, args.outdir, args.batch_size, args.binning_method, args.bigwig_paths_file, args.bins_per_gene, args.min_mean_total_covg, args.max_corr, args.max_bin_width, args.bin_width_coding, args.bin_width_noncoding, args.n_bins, args.batch)
    elif args.subcommand == 'prepare':
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

if __name__ == '__main__':
    cli()
