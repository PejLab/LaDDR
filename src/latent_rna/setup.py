"""Process annotations and determine gene batches"""

import gzip
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd

def interval_union(intervals: list[list[int]]) -> list[list[int]]:
    """Get the union of a list of intervals
    
    Args:
        intervals: List of 2-element lists

    Returns:
        List of 2-element lists
    """
    intervals = sorted(intervals)
    union = []
    for start, end in intervals:
        if not union or union[-1][1] < start:
            union.append([start, end])
        else:
            union[-1][1] = max(union[-1][1], end)
    return union

def get_exon_regions(gtf: Path) -> pd.DataFrame:
    """Get exonic regions from a GTF file

    Gets exons that are in a protein-coding gene, not in a retained intron or
    readthrough transcript, and returns the union of all exons for each gene.

    Args:
        gtf: Path to GTF file

    Returns:
        DataFrame with columns 'gene_id', 'seqname', 'strand', 'start', and
        'end'
    """
    if str(gtf).endswith('.gz'):
        opener = gzip.open(gtf, 'rt')
    else:
        opener = open(gtf, 'r')
    exclude = ['retained_intron', 'readthrough_transcript']
    gene_seqname_strand = {}
    tx_to_gene = {}
    tx_exons = defaultdict(list)
    # Record gene info, acceptable transcripts, and exons
    with opener as gtf:
        for row in gtf:
            if row[0] == '#': continue # skip header
            row = row.strip().split('\t')
            if row[2] == 'gene':
                gene_id = row[8].split('gene_id "')[1].split('"')[0]
                if ('gene_biotype "protein_coding"' in row[8] or 'gene_type "protein_coding"' in row[8]):
                        gene_seqname_strand[gene_id] = (row[0], row[6])
            elif row[2] == 'transcript':
                if any(x in row[8] for x in exclude):
                    continue
                gene_id = row[8].split('gene_id "')[1].split('"')[0]
                tx_id = row[8].split('transcript_id "')[1].split('"')[0]
                tx_to_gene[tx_id] = gene_id
            elif row[2] == 'exon':
                tx_id = row[8].split('transcript_id "')[1].split('"')[0]
                # Convert to 0-based
                start, end = int(row[3]) - 1, int(row[4])
                tx_exons[tx_id].append((start, end))
    # Get all exons for each gene
    gene_exons = defaultdict(list)
    for tx, exons in tx_exons.items():
        if tx not in tx_to_gene:
            continue
        gene_id = tx_to_gene[tx]
        gene_exons[gene_id].extend(exons)
    # Get exon regions for each gene
    exon_data = []
    for gene_id, exons in gene_exons.items():
        if gene_id not in gene_seqname_strand:
            continue
        seqname, strand = gene_seqname_strand[gene_id]
        exon_regions = interval_union(exons)
        # Add each exon region as a separate row
        for start, end in exon_regions:
            exon_data.append({
                'gene_id': gene_id,
                'seqname': seqname,
                'strand': strand,
                'start': start,
                'end': end
            })
    return pd.DataFrame(exon_data)

def gene_coordinates(exons: pd.DataFrame, chrom_lengths: dict) -> pd.DataFrame:
    """Get various useful coordinates for each gene

    Start and end are extended by 1kb on each side to include upstream and
    downstream regions.
    
    Args:
        exons: DataFrame that includes columns 'gene_id', 'seqname', 'strand',
          'start', and 'end'. Each row represents an exonic region.
        chrom_lengths: Dictionary of chromosome lengths to restrict noncoding
          extensions at the ends of chromosomes

    Returns:
        DataFrame with index 'gene_id' and columns 'seqname', 'window_start',
        'window_end', 'strand', and 'tss'
    """
    def gene_coord_single(exons: pd.DataFrame, chrom_lengths: dict) -> pd.DataFrame:
        """Get various useful coordinates for one gene"""
        seqname = exons['seqname'].iloc[0]
        strand = exons['strand'].iloc[0]
        start = exons['start'].min()
        end = exons['end'].max()
        tss = start if strand == '+' else end
        window_start = max(0, start - 1000)
        window_end = min(end + 1000, chrom_lengths[seqname])
        return pd.DataFrame({
            'seqname': seqname,
            'window_start': window_start,
            'window_end': window_end,
            'strand': strand,
            'tss': tss,
        }, index=[0])
    genes = exons.groupby('gene_id').apply(gene_coord_single, chrom_lengths, include_groups=False)
    genes = genes.reset_index(level=1, drop=True)
    genes = genes.sort_values(by='gene_id')
    return genes

def setup(gtf: Path, chrom_file: Path, batch_size: int, outdir: Path) -> pd.DataFrame:
    """Process annotations and determine gene batches

    Writes genes.tsv, exons.tsv.gz, and n_batches.txt to outdir.
    
    Args:
        gtf: Path to GTF file
        chrom_file: Path to chromosome lengths file
        batch_size: Number of genes per batch
        outdir: Path to output directory

    Returns:
        DataFrame with index 'gene_id' and columns 'seqname', 'window_start',
        'window_end', 'strand', and 'tss'
    """
    exons = get_exon_regions(gtf)
    chrom_df = pd.read_csv(
        chrom_file,
        sep='\t',
        header=None,
        usecols=[0,1], # Only use first two columns
        names=['chrom', 'length'],
        dtype={'chrom': str, 'length': int}
    )
    chrom_lengths = dict(zip(chrom_df['chrom'], chrom_df['length']))
    chroms = list(chrom_df['chrom'])

    genes = gene_coordinates(exons, chrom_lengths)
    n_batches = int(np.ceil(len(genes) / batch_size))
    genes['batch'] = [i // batch_size for i in range(len(genes))]

    outdir.mkdir(exist_ok=True)
    genes.to_csv(outdir / 'genes.tsv', sep='\t', index=True)

    exons = exons[['gene_id', 'seqname', 'start', 'end']]
    exons['seqname'] = pd.Categorical(exons['seqname'], categories=chroms, ordered=True)
    # Record chromosome order to sort bin BED files later
    exons = exons.sort_values(by=['seqname', 'gene_id'])
    exons.to_csv(outdir / 'exons.tsv.gz', sep='\t', index=False, compression='gzip')

    with open(outdir / 'n_batches.txt', 'w') as f:
        f.write(str(n_batches))

    return genes
