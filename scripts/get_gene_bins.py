"""Get gene feature bins from GTF file

Extract exons from annotation, add noncoding regions, and split into bins for
computing latent features.
"""

import argparse
from collections import defaultdict
import math
from pathlib import Path
from queue import PriorityQueue
from typing import Iterator
from bx.intervals.intersection import IntervalTree
from gtfparse import read_gtf
import numpy as np
import pandas as pd
import pyBigWig

class AdaptiveBin:
    def __init__(self, pos_start, pos_end, coverage, gene_covg, left=None):
        """
        Start and end are BED-style 0-based, e.g. for the first base of the
        genome, start=0, end=1. `gene_covg` is the sum of coverage across all
        bases in the gene for each sample, used for normalization.
        """
        self.pos_start = pos_start
        self.pos_end = pos_end
        self.coverage = coverage
        self.gene_covg = gene_covg
        self.left = left
        if left is not None:
            left.right = self
        self.right = None
        self.mean_total_covg = np.mean(self.coverage)
        self.corr_w_right = None

    def norm_coverage(self):
        return np.log2(self.coverage / self.gene_covg + 1)

    def get_corr_w_right(self):
        if self.right is not None:
            self.corr_w_right = np.corrcoef(self.norm_coverage(), self.right.norm_coverage())[0, 1]
        else:
            self.corr_w_right = 0

    def merge_with_right(self):
        self.pos_end = self.right.pos_end
        right = self.right
        self.coverage = np.sum([self.coverage, right.coverage], axis=0)
        self.right = right.right
        if self.right is not None:
            self.right.left = self
        self.mean_total_covg = np.mean(self.coverage)
        if self.corr_w_right is not None:
            self.get_corr_w_right()

    def __lt__(self, other):
        return self.mean_total_covg < other.mean_total_covg

def get_adaptive_bins(covg: np.array, min_mean_total_covg: float, max_corr: float) -> tuple:
    q = PriorityQueue()
    # For normalization, set floor to 1 to avoid division by zero
    gene_covg = np.maximum(np.sum(covg, axis=0), 1)
    for i in range(covg.shape[0]):
        left = b if i > 0 else None
        b = AdaptiveBin(i, i + 1, covg[i, :], gene_covg, left)
        if i == 0:
            start = b
        q.put((b.mean_total_covg, b))

    ## Merge low-coverage bins
    removed = set() # Keep track of bins that have been merged into others, to ignore since they can't be removed from the queue
    while True:
        lowest_covg_bin = q.get()[1]
        if lowest_covg_bin.pos_start in removed:
            continue
        if lowest_covg_bin.mean_total_covg >= min_mean_total_covg:
            break
        left, right = lowest_covg_bin.left, lowest_covg_bin.right
        if (left is not None) and (right is None or left.mean_total_covg < right.mean_total_covg):
            left.merge_with_right()
        elif right is not None:
            removed.add(right.pos_start)
            lowest_covg_bin.merge_with_right()
            q.put((lowest_covg_bin.mean_total_covg, lowest_covg_bin))
        else:
            break

    ## Normalize coverage and calculate bin correlations
    current = start
    q = PriorityQueue()
    while current is not None:
        current.get_corr_w_right()
        q.put((-current.corr_w_right, current))
        current = current.right

    while True:
        highest_corr_bin = q.get()[1]
        if highest_corr_bin.corr_w_right <= max_corr:
            break
        highest_corr_bin.merge_with_right()
        q.put((-highest_corr_bin.corr_w_right, highest_corr_bin))

    ## Assemble bins into BED format
    # df = pd.DataFrame(columns=['seqname', 'start', 'end', 'strand', 'feature'])
    starts, ends = [], []
    current = start
    while current is not None:
        # bin_id = f'{gene_info["gene_id"]}_{current.pos_start - 1000}_{current.pos_end - 1000}'
        starts.append(current.pos_start)
        ends.append(current.pos_end)
        current = current.right
    return starts, ends

def load_covg_from_bigwigs(bigwig_paths: list, seqname: str, start: int, end: int) -> np.array:
    """Load coverage data from bigWig files

    Args:
        bigwig_paths: List of paths to bigWig files to load
        seqname: Chromosome name
        start: Start position (0-based)
        end: End position (0-based)

    Returns:
        Array of shape (end - start, len(bigwig_paths)) with coverage data
    """
    covg = np.zeros((end - start, len(bigwig_paths)))
    for i, path in enumerate(bigwig_paths):
        with pyBigWig.open(str(path)) as bw:
            try:
                covg[:, i] = bw.values(seqname, start, end)
            except RuntimeError as e:
                print(f"RuntimeError for {seqname}:{start}-{end} in file {bw.stem}: {e}", flush=True)
                raise e
    return covg

def get_adaptive_bins_batch(genes: pd.DataFrame, bigwig_paths: list, min_mean_total_covg: float, max_corr: float) -> pd.DataFrame:
    """Get adaptive bins for each gene in the annotation
    
    Args:
        genes: DataFrame with index 'gene_id' and columns 'seqname', 'start',
          'end', 'strand', 'feature'
        bigwig_paths: List of paths to bigWig files to use for coverage data
        min_mean_total_covg: Minimum allowed mean total coverage per sample for
          a bin
        max_corr: Maximum allowed correlation between normalized coverage of
          adjacent bins

    Returns:
        DataFrame with columns 'gene_id', 'seqname', 'start', 'end', 'strand', 'feature'
    """
    def get_adaptive_bins_gene(gene: pd.DataFrame) -> pd.DataFrame:
        seqname, start, end, strand = gene[['seqname', 'start', 'end', 'strand']].iloc[0]
        covg = load_covg_from_bigwigs(bigwig_paths, seqname, start, end)
        starts, ends = get_adaptive_bins(covg, min_mean_total_covg, max_corr)
        bins = pd.DataFrame({
            'seqname': seqname,
            'start': starts,
            'end': ends,
            'strand': strand,
            'feature': 'adaptive',
        })
        bins['start'] = start + bins['start']
        bins['end'] = start + bins['end']
        assert bins.iloc[0].start == start 
        assert bins.iloc[-1].end == end
        assert bins.start.unique().shape[0] == bins.shape[0]
        return bins
    bins = genes.groupby('gene_id').apply(lambda x: get_adaptive_bins_gene(x), include_groups=False)
    bins = bins.reset_index()
    return bins

def add_noncoding_regions(anno: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:
    """Add noncoding regions to the annotations for one gene
    
    Inserts an intron between each pair of consecutive exons, and adds
    fixed-length (1kb) upstream and downstream regions.

    Args:
        anno: DataFrame with columns 'seqname', 'start', 'end', 'strand',
          and 'feature'. Each row represents an exonic region.
        genes: DataFrame with one row and columns 'seqname', 'start', 'end',
          'strand'. Start and end give the range to be binned, which includes
          the upstream and downstream regions.

    Returns:
        DataFrame that is the input DataFrame with additional rows for noncoding
        regions
    """
    seqname, start, end, strand = gene[['seqname', 'start', 'end', 'strand']]
    assert strand in {'+', '-'}
    assert anno.start.unique().shape[0] == anno.shape[0]
    assert (anno.feature == 'exon').all()
    features, starts, ends = [], [], []
    if anno.iloc[0].start > start:
        features.append('upstream' if strand == '+' else 'downstream')
        starts.append(start)
        ends.append(anno.iloc[0].start)
    for i in range(anno.shape[0] - 1):
        # Only add intron if there is a gap between exons:
        if anno.iloc[i + 1].start - anno.iloc[i].end > 0:
            features.append('intron')
            starts.append(anno.iloc[i].end)
            ends.append(anno.iloc[i + 1].start)
    if anno.iloc[-1].end < end:
        features.append('downstream' if strand == '+' else 'upstream')
        starts.append(anno.iloc[-1].end)
        ends.append(end)
    nc = pd.DataFrame({
        'seqname': seqname,
        'start': starts,
        'end': ends,
        'strand': strand,
        'feature': features,
    })
    anno = pd.concat([anno, nc])
    # anno = anno.sort_values(by='start')  # Sort all together for efficiency
    return anno

def split_region_bin_width(region: pd.DataFrame, bin_width_coding: int, bin_width_noncoding: int) -> pd.DataFrame:
    """Split a region into a fixed number of equal-sized bins
    
    Args:
        region: DataFrame with one row representing the region to split. Must
          have columns 'seqname', 'start', 'end', 'strand', 'feature'.
        bin_width_coding: Width of bins for coding regions (feature == 'exon')
        bin_width_noncoding: Width of bins for noncoding regions (feature !=
          'exon')

    Returns:
        DataFrame with one row per bin, with columns 'seqname', 'start', 'end',
        'strand', 'feature'
    """
    assert region.shape[0] == 1
    region = list(region.itertuples(index=False))[0]
    bin_width = bin_width_coding if region.feature == 'exon' else bin_width_noncoding
    starts = np.arange(region.start, region.end, bin_width)
    ends = starts + bin_width
    ends[-1] = region.end
    bins = pd.DataFrame({
        'seqname': region.seqname,
        'start': starts,
        'end': ends,
        'strand': region.strand,
        'feature': region.feature,
    })
    assert bins.iloc[0].start == region.start
    assert bins.iloc[-1].end == region.end
    assert bins.start.unique().shape[0] == bins.shape[0]
    return bins

def split_region_n_bins(region: pd.DataFrame, n_bins: int) -> pd.DataFrame:
    """Split a region into a fixed number of equal-sized bins

    Args:
        region: DataFrame with one row representing the region to split
        n_bins: Number of bins to split the region into

    Returns:
        DataFrame with one row per bin, with columns 'seqname', 'start', 'end',
        'strand', 'feature'
    """
    assert region.shape[0] == 1
    region = list(region.itertuples(index=False))[0]
    if region.end - region.start < n_bins:
        # If too small to split, do one bin per base and return fewer than n_bins bins
        starts = np.arange(region.start, region.end)
        bins = pd.DataFrame({
            'seqname': region.seqname,
            'start': starts,
            'end': starts + 1,
            'strand': region.strand,
            'feature': region.feature,
        })
    else:
        posns = np.linspace(region.start, region.end, n_bins + 1, dtype=int)
        bins = pd.DataFrame({
            'seqname': region.seqname,
            'start': posns[:-1],
            'end': posns[1:],
            'strand': region.strand,
            'feature': region.feature,
        })
    assert bins.iloc[0].start == region.start
    assert bins.iloc[-1].end == region.end
    assert bins.start.unique().shape[0] == bins.shape[0]
    return bins

def split_regions_bin_width(anno: pd.DataFrame, bin_width_coding: int, bin_width_noncoding: int) -> pd.DataFrame:
    """Split each region into bins of a certain width"""
    anno['start2'] = anno['start']
    anno = anno.groupby(['gene_id', 'start2'])
    bins = anno.apply(lambda x: split_region_bin_width(x, bin_width_coding, bin_width_noncoding), include_groups=False)
    bins = bins.reset_index(level='gene_id')
    bins = bins.reset_index(drop=True)
    return bins

def split_regions_n_bins(anno: pd.DataFrame, n_bins: int) -> pd.DataFrame:
    """Split each region into a fixed number of equal-sized bins"""
    anno['start2'] = anno['start']
    anno = anno.groupby(['gene_id', 'start2'])
    bins = anno.apply(lambda x: split_region_n_bins(x, n_bins), include_groups=False)
    bins = bins.reset_index(level='gene_id')
    bins = bins.reset_index(drop=True)
    return bins

def bins_overlap_exons(bins: pd.DataFrame, exons: pd.DataFrame) -> Iterator[bool]:
    """Check if each bin overlaps with an exon of another gene"""
    exon_intervals = defaultdict(IntervalTree)
    for i, exon in enumerate(exons.itertuples(index=False)):
        exon_intervals[exon.seqname].add(exon.start, exon.end, [i, exon.gene_id])
    for bin in bins.itertuples(index=False):
        overlaps = exon_intervals[bin.seqname].find(bin.start, bin.end)
        overlaps = [x for x in overlaps if x[1] != bin.gene_id]
        yield len(overlaps) > 0

def remove_bins_overlapping_exon(bins: pd.DataFrame, exons: pd.DataFrame) -> pd.DataFrame:
    """Remove bins that overlap with an exon of another gene"""
    overlap = np.array(list(bins_overlap_exons(bins, exons)))
    percentage = 100 * sum(overlap) / len(overlap)
    types_removed = bins.loc[overlap, 'feature'].value_counts()
    print(f'Removed {sum(overlap)} bins ({percentage:.6}%) that overlapped with an exon of another gene')
    print(f'Bin types removed: {types_removed.to_dict()}')
    bins = bins.loc[~overlap, :]
    return bins

def load_chromosomes(chrom_file: Path) -> pd.DataFrame:
    """Load chromosome length and sort order info"""
    chrom = pd.read_csv(
        chrom_file,
        sep='\t',
        header=None,
        names=['chrom', 'length'],
        dtype={'chrom': str, 'length': int}
    )
    return chrom

def gene_coordinates(anno: pd.DataFrame, chrom_lengths: dict) -> pd.DataFrame:
    """Get various useful coordinates for each gene

    Start and end are extended by 1kb on each side to include upstream and
    downstream regions.
    
    Args:
        anno: DataFrame that includes columns 'gene_id', 'seqname', 'start',
          'end', and 'strand'. Each row represents an exonic region.
        chrom_lengths: Dictionary of chromosome lengths to restrict noncoding
          extensions at the ends of chromosomes

    Returns:
        DataFrame with index 'gene_id' and columns 'seqname', 'start', 'end',
        'strand', and 'tss'
    """
    def gene_coord_single(anno: pd.DataFrame, chrom_lengths: dict) -> pd.DataFrame:
        """Get various useful coordinates for one gene"""
        # gene_id = anno['gene_id'].iloc[0]
        seqname = anno['seqname'].iloc[0]
        strand = anno['strand'].iloc[0]
        start = anno['start'].min()
        end = anno['end'].max()
        tss = start if strand == '+' else end
        start = max(0, start - 1000)
        end = min(end + 1000, chrom_lengths[seqname])
        return pd.DataFrame({
            'seqname': seqname,
            'start': start,
            'end': end,
            'strand': strand,
            'tss': tss,
        }, index=[0])
    genes = anno.groupby('gene_id').apply(gene_coord_single, chrom_lengths, include_groups=False)
    genes = genes.reset_index(level=1, drop=True)
    genes = genes.sort_values(by='gene_id')
    return genes

def name_bins_with_gene_coords(bins: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    """Name bins with gene-relative coordinates"""
    assert bins['gene_id'].nunique() == 1
    # assert 'exon' in bins['feature'].values
    gene_id = bins['gene_id'].iloc[0]
    tss = genes.at[gene_id, 'tss']
    # Using strand to orient correctly, define gene start as the start of the first exon
    if bins.iloc[0].strand == '+':
        starts = bins['start'] - tss
        ends = bins['end'] - tss
    else:
        starts = tss - bins['end']
        ends = tss - bins['start']
    bins['name'] = bins['gene_id'] + '_' + starts.astype(str) + '_' + ends.astype(str)
    return bins

def save_bed(bins: pd.DataFrame, chromosomes: list, genes: pd.DataFrame, outfile: Path):
    """Prepare BED format"""
    bins['gene_id2'] = bins['gene_id']
    bins = bins.groupby('gene_id2').apply(name_bins_with_gene_coords, genes, include_groups=False)
    bins = bins[['seqname', 'start', 'end', 'name', 'feature', 'strand']]
    bins['seqname'] = pd.Categorical(bins['seqname'], categories=chromosomes, ordered=True)
    bins = bins.sort_values(by=['seqname', 'start'])
    bins.to_csv(outfile, sep='\t', index=False, header=False)

parser = argparse.ArgumentParser(description='Get gene feature bins from GTF file')
parser.add_argument('-g', '--gtf', type=Path, required=True, metavar='FILE', help='Transcript annotation in GTF format. Must be the collapsed annotation produced by `collapse_annotation.py`.')
parser.add_argument('-c', '--chromosomes', type=Path, required=True, metavar='FILE', help='Chromosome lengths file, e.g. chrNameLength.txt from STAR index, to sort chromosomes.')
parser.add_argument('--outdir', type=Path, required=True, metavar='PATH', help='Directory in which to save per-batch BED files.')
parser.add_argument('--binning-method', choices=['adaptive1', 'adaptive2', 'bin-width', 'n-bins'], default='n-bins', help='Whether to determine bins adaptively from coverage data, split all coding/noncoding regions into fixed-width bins, or split each region into a fixed number of bins. (default: %(default)s)')
parser.add_argument('--bigwig-paths-file', type=Path, metavar='FILE', help='For method "adaptive", file containing list of paths to per-sample bigWig files to use for coverage data.')
parser.add_argument('--min-mean-total-covg', type=float, default=128, metavar='FLOAT', help='For method "adaptive", minimum allowed mean total coverage per sample for a bin. (default: %(default)s)')
parser.add_argument('--max-corr', type=float, default=0.8, metavar='FLOAT', help='For method "adaptive", maximum allowed correlation between normalized coverage of adjacent bins. (default: %(default)s)')
parser.add_argument('--bin-width-coding', type=int, default=16, metavar='N', help='For method "bin-width", width of bins for coding (exonic) regions of the genes. (default: %(default)s)')
parser.add_argument('--bin-width-noncoding', type=int, default=128, metavar='N', help='For method "bin-width", width of bins for noncoding (non-exonic) regions of the genes. (default: %(default)s)')
parser.add_argument('--n-bins', type=int, default=24, metavar='N', help='For method "n-bins", number of bins to split each feature (exon, intron, etc.) into. (default: %(default)s)')
parser.add_argument('--batch-size', type=int, default=200, metavar='N', help='Number of genes (at most) per batch. (default: %(default)s)')
parser.add_argument('--batch', type=int, metavar='N', help='Batch ID to process. Batch IDs are integers starting from 0. If omitted, all batches will be processed.')
args = parser.parse_args()

anno = read_gtf(args.gtf)
# Newer versions return a polars DF by default, but not all versions allow
# return type to be specified, so this handles older and newer versions:
if type(anno).__module__ == 'polars.dataframe.frame':
    anno = anno.to_pandas()

# Require only one isoform per gene:
n_iso = anno.groupby('gene_id')['transcript_id'].nunique()
assert (n_iso == 1).all(), 'Multiple isoforms per gene found in GTF file. Use `collapse_annotation.py` with the `--collapse_only` argument to collapse the annotation to a single isoform per gene.'

anno['start'] = anno['start'] - 1  # Convert to 0-based
anno = anno.loc[anno['feature'] == 'exon', :]
exons = anno.copy() # Save all exons to check for overlap later
type_col = 'gene_type' if 'gene_type' in anno.columns else 'gene_biotype'
anno = anno.loc[anno[type_col] == 'protein_coding', :]
anno = anno[['gene_id', 'seqname', 'start', 'end', 'strand', 'feature']]
anno = anno.sort_values(by=['gene_id', 'start'])

chromosomes = load_chromosomes(args.chromosomes)
chrom_lengths = dict(zip(chromosomes['chrom'], chromosomes['length']))
chromosomes = list(chromosomes['chrom'])

genes = gene_coordinates(anno, chrom_lengths)

# Split genes into batches
# genes = anno['gene_id'].sort_values().unique()
n_batches = math.ceil(genes.shape[0] / args.batch_size)
gene_batch = {}
for i, gene_id in enumerate(genes.index):
    gene_batch[gene_id] = i // args.batch_size
anno['batch'] = [gene_batch[gene_id] for gene_id in anno['gene_id']]
genes['batch'] = [gene_batch[gene_id] for gene_id in genes.index]
anno = anno.groupby('batch')
genes = genes.groupby('batch')

if args.binning_method in ['adaptive1', 'adaptive2']:
    assert args.bigwig_paths_file is not None, 'For adaptive binning, must provide --bigwig-paths-file'
    with open(args.bigwig_paths_file) as f:
        bigwig_paths = f.read().splitlines()

args.outdir.mkdir(exist_ok=True)
genes_removed = 0
batches = range(n_batches) if args.batch is None else [args.batch]
for batch_id in batches:
    batch_genes = genes.get_group(batch_id)
    if args.binning_method == 'adaptive1':
        batch_bins = get_adaptive_bins_batch(batch_genes, bigwig_paths, args.min_mean_total_covg, args.max_corr)
        batch_bins = remove_bins_overlapping_exon(batch_bins, exons)
    else:
        batch_anno = anno.get_group(batch_id)
        # For each gene, add introns, upstream, and downstream regions:
        batch_anno = batch_anno.groupby('gene_id').apply(
            lambda x: add_noncoding_regions(x, batch_genes.loc[x.name]),
            include_groups=False
        )
        batch_anno = batch_anno.reset_index(level='gene_id')
        batch_anno = batch_anno.sort_values(by=['gene_id', 'start'])
        if args.binning_method == 'bin-width':
            batch_bins = split_regions_bin_width(batch_anno, args.bin_width_coding, args.bin_width_noncoding)
        else:
            batch_bins = split_regions_n_bins(batch_anno, args.n_bins)
        batch_bins = remove_bins_overlapping_exon(batch_bins, exons)

        # Remove any genes with no exon regions remaining:
        n_genes_before = batch_bins['gene_id'].nunique()
        batch_bins = batch_bins.groupby('gene_id').filter(lambda x: 'exon' in x['feature'].values)
        genes_removed += n_genes_before - batch_bins['gene_id'].nunique()
    outfile = args.outdir / f'batch_{batch_id}.bed.gz'
    save_bed(batch_bins, chromosomes, batch_genes, outfile)
    print(f'Saved batch {batch_id} to {outfile}')

if genes_removed > 0:
    print(f'Removed {genes_removed} genes with no unique exonic regions')
