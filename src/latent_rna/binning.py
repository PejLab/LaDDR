"""Get gene feature bins from GTF file

Extract exons from annotation, add noncoding regions, and split into bins for
computing latent features.
"""

from collections import defaultdict
import math
from pathlib import Path
from queue import PriorityQueue
from typing import Iterator, Optional
from bx.intervals.intersection import IntervalTree
import numpy as np
import pandas as pd
from .coverage import base_covg_from_bigwigs

class AdaptiveBin:
    """Used for the adaptive_covgcorr binning method"""
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

def get_adaptive_bins_covgcorr(covg: np.array, min_mean_total_covg: float, max_corr: float) -> tuple:
    """Get adaptive bins for a single gene using method "adaptive_covgcorr"
    
    Args:
        covg: Coverage data for the gene, shape (n_bins, n_samples)
        min_mean_total_covg: Minimum allowed mean total coverage per sample for
          a bin
        max_corr: Maximum allowed correlation between normalized coverage of
          adjacent bins

    Returns:
        Tuple of start and end positions of bins
    """
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

    starts, ends = [], []
    current = start
    while current is not None:
        starts.append(current.pos_start)
        ends.append(current.pos_end)
        current = current.right
    return starts, ends

def get_adaptive_bins_covgcorr_batch(
        genes: pd.DataFrame,
        bigwig_paths: list[Path],
        min_mean_total_covg: float,
        max_corr: float
) -> pd.DataFrame:
    """Get adaptive bins for each gene in the annotation using method "adaptive_covgcorr"
    
    Args:
        genes: DataFrame with index 'gene_id' and columns 'seqname',
          'window_start', 'window_end', 'strand'
        bigwig_paths: List of paths to bigWig files to use for coverage data
        min_mean_total_covg: Minimum allowed mean total coverage per sample for
          a bin
        max_corr: Maximum allowed correlation between normalized coverage of
          adjacent bins

    Returns:
        DataFrame with columns 'gene_id', 'seqname', 'start', 'end', 'strand', 'feature'
    """
    def get_adaptive_bins_gene(gene: pd.DataFrame) -> pd.DataFrame:
        seqname, window_start, window_end, strand = gene[['seqname', 'window_start', 'window_end', 'strand']].iloc[0]
        covg = base_covg_from_bigwigs(bigwig_paths, seqname, window_start, window_end)
        starts, ends = get_adaptive_bins_covgcorr(covg, min_mean_total_covg, max_corr)
        bins = pd.DataFrame({
            'seqname': seqname,
            'start': starts,
            'end': ends,
            'strand': strand,
            'feature': 'adaptive',
        })
        bins['start'] = window_start + bins['start']
        bins['end'] = window_start + bins['end']
        assert bins.iloc[0].start == window_start 
        assert bins.iloc[-1].end == window_end
        assert bins.start.unique().shape[0] == bins.shape[0]
        return bins
    bins = genes.groupby('gene_id').apply(lambda x: get_adaptive_bins_gene(x), include_groups=False)
    bins = bins.reset_index()
    return bins

def estimate_var_sum_per_gene(
        genes: pd.DataFrame,
        bigwig_paths: list[Path],
        median_coverage: float,
        covg_diff: bool = False,
        pseudocount: float = 8
) -> pd.DataFrame:
    """Use a subsample of genes to estimate mean sum of variance across a gene

    This uses either the variance of log-coverage, or the variance of the diff
    in log-coverage between adjacent positions. By dividing this by the desired
    number of bins per gene, we can estimate the cumulative variance threshold
    per bin that produces the desired number of bins on average.
    
    Args:
        genes: DataFrame with index 'gene_id' and columns 'seqname',
          'window_start', 'window_end', 'strand',
        bigwig_paths: List of paths to bigWig files to use for coverage data
        covg_diff: If True, use the difference in coverage between adjacent bins
        pseudocount: Pseudocount to add to coverage values (mean coverage per
          bp for each bin)

    Returns:
        Mean sum of variance per bin across all genes
    """
    total = 0
    for gene in genes.itertuples(index=False):
        seqname, window_start, window_end = gene.seqname, gene.window_start, gene.window_end
        covg = base_covg_from_bigwigs(bigwig_paths, seqname, window_start, window_end, median_coverage)
        covg = np.log2(covg + pseudocount)
        if covg_diff:
            covg = np.diff(covg, axis=0, append=np.log2(pseudocount))
        total += np.sum(np.var(covg, axis=1))
    return total / genes.shape[0]

def variance_threshold(
        genes: pd.DataFrame,
        bigwig_paths: list[Path],
        bins_per_gene: int,
        median_coverage: float = None,
        covg_diff: bool = False,
        pseudocount: float = 8
) -> float:
    """Estimate the variance threshold needed per bin

    Due to spikes in variance, especially for diff of log-coverage, a basic
    threshold calculation based on total cumulative variance per gene and
    desired number of bins per gene will produce too few bins. This function
    calculates an initial threshold, calculates bins on the sample of genes, and
    then calibrates the threshold based on the actual number of bins produced.

    Args:
        genes: DataFrame with index 'gene_id' and columns 'seqname',
          'window_start', 'window_end', 'strand',
        bigwig_paths: List of paths to bigWig files to use for coverage data
        bins_per_gene: Desired number of bins per gene on average
        median_coverage: Median of sumData across all samples. If provided,
          coverage values will be scaled by sumData/median_coverage to normalize
          for sequencing depth. If None, no scaling is applied.
        covg_diff: If True, use the difference in coverage between adjacent bins
        pseudocount: Pseudocount to add to coverage values (mean coverage per
          bp for each bin)

    Returns:
        Calibrated variance threshold per bin
    """
    n_genes = min(128, genes.shape[0])
    sub_genes = genes.sample(n_genes)
    var_per_gene = estimate_var_sum_per_gene(sub_genes, bigwig_paths, median_coverage, covg_diff, pseudocount)
    initial_threshold = var_per_gene / bins_per_gene
    print(f'Initial threshold: {initial_threshold}')
    bins = get_adaptive_bins_var_batch(sub_genes, bigwig_paths, initial_threshold, covg_diff, pseudocount)
    actual_bins_per_gene = bins.shape[0] / n_genes
    new_threshold = initial_threshold * actual_bins_per_gene / bins_per_gene
    print(f'Adjusted threshold: {new_threshold}')
    return new_threshold

def get_adaptive_bins_var_batch(
        genes: pd.DataFrame,
        bigwig_paths: list[Path],
        var_per_bin: float,
        covg_diff: bool = False,
        pseudocount: float = 8
) -> pd.DataFrame:
    """Get adaptive bins for all genes based on coverage variance
    
    Args:
        genes: DataFrame with index 'gene_id' and columns 'seqname',
          'window_start', 'window_end', 'strand',
        bigwig_paths: List of paths to bigWig files to use for coverage data
        var_per_bin: Approximate sum of per-bp variance of log-coverage across
          samples for each bin
        covg_diff: If True, use the difference in coverage between adjacent bins
        pseudocount: Pseudocount to add to coverage values (mean coverage per
          bp for each bin)

    Returns:
        DataFrame with columns 'gene_id', 'seqname', 'start', 'end', 'strand', 'feature'
    """
    def get_adaptive_bins_var(covg: np.array, var_per_bin: float) -> tuple:
        covg = np.log2(covg + pseudocount)
        if covg_diff:
            covg = np.diff(covg, axis=0, append=np.log2(pseudocount))
        variance = np.var(covg, axis=1)
        starts, ends = [0], []
        varsum = 0
        for i, var in enumerate(variance[:-1]):
            varsum += var
            if varsum >= var_per_bin:
                ends.append(i + 1)
                starts.append(i + 1)
                varsum = 0
        ends.append(covg.shape[0])
        return starts, ends
    def get_adaptive_bins_var_gene(gene: pd.DataFrame) -> pd.DataFrame:
        seqname, window_start, window_end, strand = gene[['seqname', 'window_start', 'window_end', 'strand']].iloc[0]
        covg = base_covg_from_bigwigs(bigwig_paths, seqname, window_start, window_end)
        starts, ends = get_adaptive_bins_var(covg, var_per_bin)
        bins = pd.DataFrame({
            'seqname': seqname,
            'start': starts,
            'end': ends,
            'strand': strand,
            'feature': 'adaptive',
        })
        bins['start'] = window_start + bins['start']
        bins['end'] = window_start + bins['end']
        assert bins.iloc[0].start == window_start 
        assert bins.iloc[-1].end == window_end
        assert bins.start.unique().shape[0] == bins.shape[0]
        return bins
    bins = genes.groupby('gene_id').apply(lambda x: get_adaptive_bins_var_gene(x), include_groups=False)
    bins = bins.reset_index()
    return bins

def add_noncoding_regions(anno: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:
    """Add noncoding regions to the annotations for one gene
    
    Inserts an intron between each pair of consecutive exons, and adds up to 1kb
    upstream and downstream regions.

    Args:
        anno: DataFrame with columns 'seqname', 'start', 'end', 'strand',
          and 'feature'. Each row represents an exonic region.
        gene: DataFrame with one row and columns 'seqname', 'window_start',
          'window_end', 'strand'. Start and end give the range to be binned, which includes
          the upstream and downstream regions.

    Returns:
        DataFrame that is the input DataFrame with additional rows for noncoding
        regions
    """
    seqname, window_start, window_end, strand = gene[['seqname', 'window_start', 'window_end', 'strand']]
    assert strand in {'+', '-'}
    assert anno.start.unique().shape[0] == anno.shape[0]
    assert (anno.feature == 'exon').all()
    features, starts, ends = [], [], []
    if anno.iloc[0].start > window_start:
        features.append('upstream' if strand == '+' else 'downstream')
        starts.append(window_start)
        ends.append(anno.iloc[0].start)
    for i in range(anno.shape[0] - 1):
        # Only add intron if there is a gap between exons:
        if anno.iloc[i + 1].start - anno.iloc[i].end > 0:
            features.append('intron')
            starts.append(anno.iloc[i].end)
            ends.append(anno.iloc[i + 1].start)
    if anno.iloc[-1].end < window_end:
        features.append('downstream' if strand == '+' else 'upstream')
        starts.append(anno.iloc[-1].end)
        ends.append(window_end)
    nc = pd.DataFrame({
        'seqname': seqname,
        'start': starts,
        'end': ends,
        'strand': strand,
        'feature': features,
    })
    anno = pd.concat([anno, nc])
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
    """Split each region into bins of a certain width
    
    Args:
        anno: DataFrame with columns 'gene_id', 'seqname', 'start', 'end',
          'strand', and 'feature'.
        bin_width_coding: Width of bins for coding regions (feature == 'exon')
        bin_width_noncoding: Width of bins for noncoding regions (feature !=
          'exon')

    Returns:
        DataFrame with one row per bin, with columns 'seqname', 'start', 'end',
        'strand', 'feature'
    """
    anno['start2'] = anno['start']
    anno = anno.groupby(['gene_id', 'start2'])
    bins = anno.apply(lambda x: split_region_bin_width(x, bin_width_coding, bin_width_noncoding), include_groups=False)
    bins = bins.reset_index(level='gene_id')
    bins = bins.reset_index(drop=True)
    return bins

def split_regions_n_bins(anno: pd.DataFrame, n_bins: int) -> pd.DataFrame:
    """Split each region into a fixed number of equal-sized bins
    
    Args:
        anno: DataFrame with columns 'gene_id', 'seqname', 'start', 'end',
          'strand', and 'feature'.
        n_bins: Number of bins to split each region into

    Returns:
        DataFrame with one row per bin, with columns 'seqname', 'start', 'end',
        'strand', 'feature'
    """
    anno['start2'] = anno['start']
    anno = anno.groupby(['gene_id', 'start2'])
    bins = anno.apply(lambda x: split_region_n_bins(x, n_bins), include_groups=False)
    bins = bins.reset_index(level='gene_id')
    bins = bins.reset_index(drop=True)
    return bins

def split_large_bins(anno: pd.DataFrame, max_bin_width: int) -> pd.DataFrame:
    """Split large bins into smaller bins
    
    Args:
        anno: DataFrame representing the bins to split. Must have columns
          'seqname', 'start', 'end', 'strand', 'feature'.
        max_bin_width: Maximum allowed width for a bin

    Returns:
        DataFrame with one row per bin, with columns 'seqname', 'start', 'end',
        'strand', 'feature'
    """
    def split_large_bin(bin: pd.DataFrame) -> pd.DataFrame:
        assert bin.shape[0] == 1
        bin2 = list(bin.itertuples(index=False))[0]
        if bin2.end - bin2.start <= max_bin_width:
            return bin
        bin = bin2
        n_bins = math.ceil((bin.end - bin.start) / max_bin_width)
        posns = np.linspace(bin.start, bin.end, n_bins + 1, dtype=int)
        bins = pd.DataFrame({
            'seqname': bin.seqname,
            'start': posns[:-1],
            'end': posns[1:],
            'strand': bin.strand,
            'feature': bin.feature,
        })
        assert bins.iloc[0].start == bin.start
        assert bins.iloc[-1].end == bin.end
        assert bins.start.unique().shape[0] == bins.shape[0]
        return bins
    anno['start2'] = anno['start']
    anno = anno.groupby(['gene_id', 'start2'])
    bins = anno.apply(lambda x: split_large_bin(x), include_groups=False)
    bins = bins.reset_index(level='gene_id')
    bins = bins.reset_index(drop=True)
    return bins


def bins_overlap_exons(bins: pd.DataFrame, exons: pd.DataFrame) -> Iterator[bool]:
    """Check if each bin overlaps with an exon of another gene
    
    Args:
        bins: DataFrame with columns 'seqname', 'start', 'end', 'strand',
          'feature', and 'gene_id'.
        exons: DataFrame with columns 'gene_id', 'seqname', 'start', and 'end'.

    Returns:
        Iterator of booleans indicating whether each bin overlaps with an exon
        of another gene
    """
    exon_intervals = defaultdict(IntervalTree)
    for i, exon in enumerate(exons.itertuples(index=False)):
        exon_intervals[exon.seqname].add(exon.start, exon.end, [i, exon.gene_id])
    for bin in bins.itertuples(index=False):
        overlaps = exon_intervals[bin.seqname].find(bin.start, bin.end)
        overlaps = [x for x in overlaps if x[1] != bin.gene_id]
        yield len(overlaps) > 0

def remove_bins_overlapping_exon(bins: pd.DataFrame, exons: pd.DataFrame) -> pd.DataFrame:
    """Remove bins that overlap with an exon of another gene
    
    Args:
        bins: DataFrame with columns 'seqname', 'start', 'end', 'strand',
          'feature', and 'gene_id'.
        exons: DataFrame with columns 'gene_id', 'seqname', 'start', and 'end'.

    Returns:
        DataFrame of filtered bins
    """
    overlap = np.array(list(bins_overlap_exons(bins, exons)))
    percentage = 100 * sum(overlap) / len(overlap)
    # types_removed = bins.loc[overlap, 'feature'].value_counts()
    print(f'Removed {sum(overlap)} bins ({percentage:.6}%) that overlapped with an exon of another gene')
    # print(f'Bin types removed: {types_removed.to_dict()}')
    bins = bins.loc[~overlap, :]
    return bins

def name_bins_with_gene_coords(bins: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    """Name bins with gene-relative coordinates
    
    Args:
        bins: DataFrame with columns 'seqname', 'start', 'end', 'strand',
          'feature', and 'gene_id'.
        genes: DataFrame with columns 'gene_id', 'seqname', 'window_start',
          'window_end', 'strand', and 'batch'.

    Returns:
        DataFrame with one row per bin, with columns 'seqname', 'start', 'end',
        'strand', 'feature', 'gene_id', and 'name'.
    """
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
    """Prepare BED format for a batch of bins
    
    Args:
        bins: DataFrame with columns 'seqname', 'start', 'end', 'strand',
          'feature', and 'gene_id'.
        chromosomes: List of chromosome names
        genes: DataFrame with columns 'gene_id', 'seqname', 'window_start',
          'window_end', 'strand', and 'batch'.
        outfile: Path to save the BED file
    """
    bins['gene_id2'] = bins['gene_id']
    bins = bins.groupby('gene_id2').apply(name_bins_with_gene_coords, genes, include_groups=False)
    bins = bins[['seqname', 'start', 'end', 'name', 'feature', 'strand']]
    bins['seqname'] = pd.Categorical(bins['seqname'], categories=chromosomes, ordered=True)
    bins = bins.sort_values(by=['seqname', 'start'])
    bins.to_csv(outfile, sep='\t', index=False, header=False)

def adaptive_binning(
    gene_file: Path,
    exon_file: Path,
    outdir: Path,
    bigwig_paths: list[Path],
    batch: Optional[int] = None,
    binning_method: str = 'adaptive_diffvar',
    var_threshold: Optional[int] = None,
    min_mean_total_covg: Optional[float] = None,
    max_corr: Optional[float] = None,
    max_bin_width: Optional[int] = None,
):
    """Adaptive binning based on RNA-seq coverage patterns

    Writes a BED file to `outdir`/`batch_id`.bed.gz for one or more batches.
    
    Args:
        gene_file: Path to a tab-delimited file with columns 'gene_id',
          'seqname', 'window_start', 'window_end', 'strand', and 'batch'.
        exon_file: Path to a tab-delimited file with columns 'gene_id',
          'seqname', 'start', and 'end'. Additional columns are ignored.
        outdir: Directory to save output BED files
        batch: Batch number to run, or None to run all batches
        bigwig_paths: List of paths to bigWig files
        binning_method: Either 'adaptive_covgcorr', 'adaptive_covgvar', or
          'adaptive_diffvar'
        var_threshold: For 'adaptive_covgvar' or 'adaptive_diffvar', the
          variance threshold for determining bin boundaries, computed to produce
          an approximate desired number of bins per gene.
        min_mean_total_covg: For 'adaptive_covgcorr', minimum allowed mean total
          coverage per sample for a bin.
        max_corr: For 'adaptive_covgcorr', maximum allowed correlation between
          normalized coverage of adjacent bins.
        max_bin_width: Maximum allowed width for a bin
    """
    genes = pd.read_csv(gene_file, sep='\t', index_col=0, dtype={'seqname': str})
    exons = pd.read_csv(exon_file, sep='\t', dtype={'seqname': str})
    chromosomes = list(exons['seqname'].unique())

    if binning_method in {'adaptive_covgvar', 'adaptive_diffvar'}:
        assert var_threshold is not None

    outdir.mkdir(exist_ok=True)
    batches = sorted(genes['batch'].unique()) if batch is None else [batch]
    for batch_id in batches:
        batch_genes = genes.groupby('batch').get_group(batch_id)
        if binning_method == 'adaptive_covgcorr':
            batch_bins = get_adaptive_bins_covgcorr_batch(batch_genes, bigwig_paths, min_mean_total_covg, max_corr)
        elif binning_method == 'adaptive_covgvar':
            batch_bins = get_adaptive_bins_var_batch(batch_genes, bigwig_paths, var_threshold)
        elif binning_method == 'adaptive_diffvar':
            batch_bins = get_adaptive_bins_var_batch(batch_genes, bigwig_paths, var_threshold, covg_diff=True)
        else:
            raise ValueError(f'Invalid binning method: {binning_method}. Expected one of: adaptive_covgcorr, adaptive_covgvar, adaptive_diffvar')
        batch_bins = split_large_bins(batch_bins, max_bin_width)
        batch_bins = remove_bins_overlapping_exon(batch_bins, exons)
        outfile = outdir / f'batch_{batch_id}.bed.gz'
        save_bed(batch_bins, chromosomes, batch_genes, outfile)
        print(f'Saved batch {batch_id} to {outfile}')

def fixed_binning(
    gene_file: Path,
    exon_file: Path,
    outdir: Path,
    batch: Optional[int] = None,
    binning_method: Optional[str] = 'fixed_width',
    max_bin_width: Optional[int] = None,
    bin_width_coding: Optional[int] = None,
    bin_width_noncoding: Optional[int] = None,
    n_bins_per_region: Optional[int] = None,
):
    """Binning based on annotated regions

    Writes a BED file to `outdir`/`batch_id`.bed.gz for one or more batches.
    
    Args:
        gene_file: Path to a tab-delimited file with columns 'gene_id', 'seqname',
          'window_start', 'window_end', 'strand', and 'batch'.
        exon_file: Path to a tab-delimited file with columns 'gene_id', 'seqname',
          'start', and 'end'. Additional columns are ignored.
        outdir: Directory to save output BED files
        batch: Batch number to run, or None to run all batches
        binning_method: Either 'fixed_width' or 'fixed_count'
        max_bin_width: Maximum allowed width for a bin
        bin_width_coding: For 'fixed_width', width of bins for coding (exonic)
          regions of the genes.
        bin_width_noncoding: For 'fixed_width', width of bins for noncoding
          (non-exonic) regions of the genes.
        n_bins_per_region: For 'fixed_count', number of bins to create per
          region.
    """
    genes = pd.read_csv(gene_file, sep='\t', index_col=0, dtype={'seqname': str})
    exons = pd.read_csv(exon_file, sep='\t', dtype={'seqname': str})
    chromosomes = list(exons['seqname'].unique())

    exons['feature'] = 'exon'
    exons['batch'] = exons['gene_id'].map(genes['batch'])
    exons_grp = exons.groupby('batch')

    outdir.mkdir(exist_ok=True)
    genes_removed = 0
    batches = range(genes['batch'].nunique()) if batch is None else [batch]
    for batch_id in batches:
        batch_genes = genes.groupby('batch').get_group(batch_id)
        batch_exons = exons_grp.get_group(batch_id)
        # For each gene, add introns, upstream, and downstream regions:
        batch_exons = batch_exons.groupby('gene_id').apply(
            lambda x: add_noncoding_regions(x, batch_genes.loc[x.name]),
            include_groups=False
        )
        batch_exons = batch_exons.reset_index(level='gene_id')
        batch_exons = batch_exons.sort_values(by=['gene_id', 'start'])
        if binning_method == 'fixed_width':
            batch_bins = split_regions_bin_width(batch_exons, bin_width_coding, bin_width_noncoding)
        else:
            batch_bins = split_regions_n_bins(batch_exons, n_bins_per_region)
        batch_bins = split_large_bins(batch_bins, max_bin_width)
        batch_bins = remove_bins_overlapping_exon(batch_bins, exons)

        # Remove any genes with no exon regions remaining:
        n_genes_before = batch_bins['gene_id'].nunique()
        batch_bins = batch_bins.groupby('gene_id').filter(lambda x: 'exon' in x['feature'].values)
        genes_removed += n_genes_before - batch_bins['gene_id'].nunique()
        outfile = outdir / f'batch_{batch_id}.bed.gz'
        save_bed(batch_bins, chromosomes, batch_genes, outfile)
        print(f'Saved batch {batch_id} to {outfile}')
    if genes_removed > 0:
        print(f'Removed {genes_removed} genes with no unique exonic regions')
