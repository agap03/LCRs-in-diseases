"""Convert protein positions of MANE Select transcripts to GRCh38 coordinates.

Helper module used by check_mane.py (and later analysis steps). It works on the
CDS segments written by prepare_mane.py (mane_cds_segments.tsv).

Coordinate conventions used here:
  - GTF / segments file: 1-based, inclusive (start..end), as in the GTF.
  - "BED intervals":     0-based, half-open [start, end), as in BED files and
                         as expected by pyfaidx slicing genome[chrom][start:end].
  - Protein positions:   1-based, inclusive (residue 1 = first amino acid).
  - "Coding order":      the order in which the segments are read by the
                         ribosome, i.e. 5'->3' along the transcript. On the minus
                         strand this is the reverse of the genomic order.
"""
import pandas as pd
from Bio.Seq import Seq


def load_segments(path="mane_cds_segments.tsv"):
    """{transcript: [(chrom, start, end, strand), ...]} in coding order (5'->3').

    Coordinates stay as in the GTF (1-based, inclusive).
    """
    seg = pd.read_csv(path, sep="\t")
    out = {}
    for tx, g in seg.groupby("transcript_id"):
        # The file is sorted by genomic start (ascending), so for plus-strand
        # transcripts this is already coding order.
        rows = list(g[["chrom", "start", "end", "strand"]].itertuples(index=False, name=None))
        # Minus-strand transcripts are read from the highest coordinate down,
        # so reverse the list to get the first coding segment first.
        out[tx] = rows[::-1] if rows[0][3] == "-" else rows
    return out


def protein_to_genome(segments, tx, aa_start, aa_end):
    """Residues aa_start..aa_end (1-based, inclusive) -> BED intervals in coding order.

    A protein region can span several exons, so the result is a list of
    intervals, one per CDS segment it touches.
    """
    # Residue k is encoded by CDS nucleotides 3k-2 .. 3k (1-based along the CDS).
    # n1..n2 is therefore the nucleotide range covering the whole region.
    n1, n2 = (aa_start - 1) * 3 + 1, aa_end * 3
    out = []
    offset = 0      # number of coding nucleotides in all segments before this one
    for chrom, start, end, strand in segments[tx]:
        length = end - start + 1
        # This segment holds CDS nucleotides offset+1 .. offset+length.
        # Intersect that with the requested range n1..n2.
        lo, hi = max(n1, offset + 1), min(n2, offset + length)
        if lo <= hi:   # the region overlaps this segment
            if strand == "+":
                # Plus strand: CDS position grows with genomic position,
                # counting from the segment start.
                g1, g2 = start + lo - offset - 1, start + hi - offset - 1
            else:
                # Minus strand: CDS position grows as genomic position falls,
                # counting from the segment end; lo maps to the higher coordinate.
                g1, g2 = end - (hi - offset - 1), end - (lo - offset - 1)
            # g1..g2 is 1-based inclusive; convert to 0-based half-open (BED).
            out.append((chrom, g1 - 1, g2, strand))
        offset += length
        if offset >= n2:   # all requested nucleotides already covered
            break
    return out


def coding_dna(genome, intervals):
    """Coding-strand DNA of BED intervals in coding order (genome = pyfaidx.Fasta).
    Minus-strand pieces are reverse-complemented one by one, then joined."""
    pieces = []
    for chrom, start, end, strand in intervals:
        # pyfaidx uses 0-based half-open slicing, matching BED intervals.
        # upper() removes soft-masking (lower-case repeats in hg38.fa).
        seq = str(genome[chrom][start:end]).upper()
        # The genome FASTA is the plus strand; for minus-strand genes take the
        # reverse complement to get the sequence as it appears in the mRNA.
        pieces.append(str(Seq(seq).reverse_complement()) if strand == "-" else seq)
    # Intervals are already in coding order, so simple concatenation gives
    # the coding sequence 5'->3'.
    return "".join(pieces)
