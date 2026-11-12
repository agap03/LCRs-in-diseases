#!/usr/bin/env python3
"""
Build a CDS segment table for MANE Select transcripts and provide a function
that converts protein residue ranges into genomic intervals (BED, 0-based).
"""
import pandas as pd

GTF = "mane_primary.gtf"
SELECT = "select_enst.txt"
OUT = "mane_cds_segments.tsv"


def build_segments():
    select = set(open(SELECT).read().split())
    gtf = pd.read_csv(GTF, sep="\t", comment="#", header=None,
                      names=["chrom", "src", "feature", "start", "end",
                             "score", "strand", "frame", "attr"])
    cds = gtf[gtf["feature"] == "CDS"].copy()
    cds["transcript_id"] = cds["attr"].str.extract(r'transcript_id "([^"]+)"')
    cds = cds[cds["transcript_id"].isin(select)]
    cds = cds[["transcript_id", "chrom", "start", "end", "strand"]]
    cds = cds.sort_values(["transcript_id", "start"]).reset_index(drop=True)
    cds.to_csv(OUT, sep="\t", index=False)
    return cds


def load_segments():
    """Return {transcript_id: list of (chrom, start, end, strand) in coding order}."""
    seg = pd.read_csv(OUT, sep="\t")
    out = {}
    for tx, g in seg.groupby("transcript_id"):
        rows = list(g[["chrom", "start", "end", "strand"]].itertuples(index=False, name=None))
        if rows[0][3] == "-":
            rows = rows[::-1]            # minus strand: coding order runs right to left
        out[tx] = rows
    return out


def protein_to_genome(segments, tx, aa_start, aa_end):
    """
    Convert residues aa_start..aa_end (1-based, inclusive) of transcript tx
    into a list of genomic intervals (chrom, bed_start, bed_end, strand).
    """
    n1, n2 = (aa_start - 1) * 3 + 1, aa_end * 3     # CDS nucleotide range, 1-based
    out, offset = [], 0
    for chrom, start, end, strand in segments[tx]:
        length = end - start + 1
        lo, hi = max(n1, offset + 1), min(n2, offset + length)
        if lo <= hi:
            if strand == "+":
                g1, g2 = start + (lo - offset - 1), start + (hi - offset - 1)
            else:
                g1, g2 = end - (hi - offset - 1), end - (lo - offset - 1)
            out.append((chrom, g1 - 1, g2, strand))      # BED: 0-based start
        offset += length
        if offset >= n2:
            break
    return out


if __name__ == "__main__":
    cds = build_segments()
    print(cds["transcript_id"].nunique(), "transcripts,", len(cds), "CDS segments ->", OUT)


"""19229 transcripts, 191809 CDS segments -> mane_cds_segments.tsv"""