"""Get gene feature bins from GTF file

Extract exons from annotation, add noncoding regions, and split into bins for
computing latent features.
"""

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Iterator
from bx.intervals.intersection import IntervalTree
from gtfparse import read_gtf
import numpy as np
import pandas as pd

def add_noncoding_regions(anno: pd.DataFrame) -> pd.DataFrame:
    """Add noncoding regions to the annotations for one gene
    
    Inserts an intron between each pair of consecutive exons, and adds
    fixed-length upstream and downstream regions.
    """
    strand = anno.iloc[0].strand
    assert strand in {'+', '-'}
    assert anno.start.unique().shape[0] == anno.shape[0]
    features, starts, ends = [], [], []
    features.append('upstream' if strand == '+' else 'downstream')
    starts.append(max(1, anno.iloc[0].start - 1000))
    ends.append(anno.iloc[0].start - 1)
    for i in range(anno.shape[0] - 1):
        # Only add intron if there is a gap between exons:
        if anno.iloc[i + 1].start - anno.iloc[i].end > 1:
            features.append('intron')
            starts.append(anno.iloc[i].end + 1)
            ends.append(anno.iloc[i + 1].start - 1)
    features.append('downstream' if strand == '+' else 'upstream')
    starts.append(anno.iloc[-1].end + 1)
    ends.append(anno.iloc[-1].end + 1000)
    nc = pd.DataFrame({
        'seqname': anno.iloc[0].seqname,
        'start': starts,
        'end': ends,
        'strand': strand,
        'feature': features,
    })
    anno = pd.concat([anno, nc])
    # anno = anno.sort_values(by='start')  # Sort all together for efficiency
    return anno

def split_region_bin_width(region: pd.DataFrame, bin_width_coding: int, bin_width_noncoding: int) -> pd.DataFrame:
    """Split a region into a fixed number of equal-sized bins"""
    assert region.shape[0] == 1
    region = list(region.itertuples(index=False))[0]
    bin_width = bin_width_coding if region.feature == 'exon' else bin_width_noncoding
    starts = np.arange(region.start, region.end + 1, bin_width)
    ends = starts + bin_width - 1
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
    """Split a region into a fixed number of equal-sized bins"""
    assert region.shape[0] == 1
    region = list(region.itertuples(index=False))[0]
    if region.end - region.start < n_bins:
        # If too small to split, do one bin per base and return fewer than n_bins bins
        starts = np.arange(region.start, region.end + 1)
        bins = pd.DataFrame({
            'seqname': region.seqname,
            'start': starts,
            'end': starts,
            'strand': region.strand,
            'feature': region.feature,
        })
    else:
        posns = np.linspace(region.start, region.end + 1, n_bins + 1, dtype=int)
        bins = pd.DataFrame({
            'seqname': region.seqname,
            'start': posns[:-1],
            'end': posns[1:] - 1,
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
        exon_intervals[exon.seqname].add(exon.start, exon.end + 1, [i, exon.gene_id]) # Interval excludes end, i.e. [start,end)
    for bin in bins.itertuples(index=False):
        overlaps = exon_intervals[bin.seqname].find(bin.start, bin.end + 1)
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

def name_bins_with_gene_coords(bins: pd.DataFrame) -> pd.DataFrame:
    """Name bins with gene-relative coordinates"""
    assert bins['gene_id'].nunique() == 1
    assert 'exon' in bins['feature'].values
    # Using strand to orient correctly, define gene start as the start of the first exon
    if bins.iloc[0].strand == '+':
        gene_start = bins.loc[bins['feature'] == 'exon', 'start'].min()
        starts = bins['start'] - gene_start
        ends = bins['end'] - gene_start
    else:
        gene_start = bins.loc[bins['feature'] == 'exon', 'end'].max()
        starts = gene_start - bins['end']
        ends = gene_start - bins['start']
    bins['name'] = bins['gene_id'] + '_' + starts.astype(str) + '_' + ends.astype(str)
    return bins

parser = argparse.ArgumentParser(description='Get gene feature bins from GTF file')
parser.add_argument('-g', '--gtf', type=Path, required=True, metavar='FILE', help='Transcript annotation in GTF format. Must be the collapsed annotation produced by `collapse_annotation.py`.')
parser.add_argument('-c', '--chromosomes', type=Path, required=True, metavar='FILE', help='Chromosome lengths file, e.g. chrNameLength.txt from STAR index, to sort chromosomes.')
parser.add_argument('-o', '--output', type=Path, required=True, metavar='FILE', help='Output file (BED).')
parser.add_argument('--binning-method', choices=['bin-width', 'n-bins'], default='bin-width', help='Whether to split all coding/noncoding regions into fixed-width bins, or split each region into a fixed number of bins. (default: %(default)s)')
parser.add_argument('--bin-width-coding', type=int, default=8, metavar='N', help='For method "bin-width", width of bins for coding (exonic) regions of the genes. (default: %(default)s)')
parser.add_argument('--bin-width-noncoding', type=int, default=16, metavar='N', help='For method "bin-width", width of bins for noncoding (non-exonic) regions of the genes. (default: %(default)s)')
parser.add_argument('--n-bins', type=int, default=10, metavar='N', help='For method "n-bins", number of bins to split each feature (exon, intron, etc.) into. (default: %(default)s)')
args = parser.parse_args()

anno = read_gtf(args.gtf)
# Newer versions return a polars DF by default, but not all versions allow
# return type to be specified, so this handles older and newer versions:
if type(anno).__module__ == 'polars.dataframe.frame':
    anno = anno.to_pandas()

# Require only one isoform per gene:
n_iso = anno.groupby('gene_id')['transcript_id'].nunique()
assert (n_iso == 1).all(), 'Multiple isoforms per gene found in GTF file. Use `collapse_annotation.py` with the `--collapse_only` argument to collapse the annotation to a single isoform per gene.'

anno = anno.loc[anno['feature'] == 'exon', :]
exons = anno.copy() # Save all exons to check for overlap later

# Keep only genes with biotype 'protein_coding':
anno = anno.loc[anno['gene_type'] == 'protein_coding', :]

anno = anno[['gene_id', 'seqname', 'start', 'end', 'strand', 'feature']]
anno = anno.sort_values(by=['gene_id', 'start'])

# For each gene, add introns, upstream, and downstream regions:
anno = anno.groupby('gene_id').apply(add_noncoding_regions, include_groups=False)
anno = anno.reset_index(level='gene_id')
anno = anno.sort_values(by=['gene_id', 'start'])
if args.binning_method == 'bin-width':
    bins = split_regions_bin_width(anno, args.bin_width_coding, args.bin_width_noncoding)
else:
    bins = split_regions_n_bins(anno, args.n_bins)
bins = remove_bins_overlapping_exon(bins, exons)

# Remove any genes with no exon regions remaining:
n_genes_before = bins['gene_id'].nunique()
bins = bins.groupby('gene_id').filter(lambda x: 'exon' in x['feature'].values)
n_genes_after = bins['gene_id'].nunique()
if n_genes_before != n_genes_after:
    print(f'Removed {n_genes_before - n_genes_after} genes with no unique exonic regions')

# Prepare BED format:
bins['start'] = bins['start'] - 1
bins['gene_id2'] = bins['gene_id']
bins = bins.groupby('gene_id2').apply(name_bins_with_gene_coords, include_groups=False)
bins = bins[['seqname', 'start', 'end', 'name', 'feature', 'strand']]

# Sort bins using chromosome order from the chromosome lengths file:
chrom = pd.read_csv(
    args.chromosomes,
    sep='\t',
    header=None,
    names=['chrom', 'length'],
    dtype={'chrom': str, 'length': int}
)
chrom = list(chrom['chrom'])
bins['seqname'] = pd.Categorical(bins['seqname'], categories=chrom, ordered=True)
bins = bins.sort_values(by=['seqname', 'start'])

bins.to_csv(args.output, sep='\t', index=False, header=False)
