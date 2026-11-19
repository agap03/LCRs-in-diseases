"""
Turn the PlaToLoCo results into one LCR table and genomic BED files.
cds_mapping.py must be in the same folder as this script.

Usage:
    python3 build_lcr_table.py full_run ~/data/mane
    python3 build_lcr_table.py full_run ~/data/mane --out-dir results --methods CAST fLPS_strict
    python3 build_lcr_table.py -h        # list all parameters

What it does, per method (default: SEG_strict, SEG_intermediate, CAST, fLPS_strict):
  - merges overlapping regions within a protein (fLPS in particular reports
    overlapping single- and multi-residue biases for the same stretch), so each
    residue is either in an LCR or not
  - describes each merged region: length, sequence, dominant residue, purity
    (fraction of the dominant residue), the residues making up >=20 %, and
    Shannon entropy
  - converts each region to GRCh38 coordinates, and makes the complementary
    "non-LCR" part of the same proteins, the comparison region for later tests

Outputs (in --out-dir, default the current folder):
  lcr_regions.tsv          one row per merged LCR per method, protein and
                           genomic description
  bed/lcr_<method>.bed     LCR genomic intervals (name = lcr_id)
  bed/nonlcr_<method>.bed  rest of the coding sequence of all analysed proteins
  lcr_summary.txt          counts per method and positive-control check
"""
import argparse
import math
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cds_mapping import load_segments, protein_to_genome  # noqa: E402

DEFAULT_METHODS = ["SEG_strict", "SEG_intermediate", "CAST", "fLPS_strict"]

# Known disease-associated repeats that the LCR detection should find:
# gene -> amino acid of the repeat (Q = polyglutamine, A = polyalanine)
CONTROLS = {"HTT": "Q", "AR": "Q", "ATXN1": "Q", "ATXN3": "Q", "TBP": "Q", "PABPN1": "A"}


def parse_args():
    p = argparse.ArgumentParser(
        description="Build the LCR table and genomic BED files from PlaToLoCo results.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("results_dir", type=Path,
                   help="folder with the merged PlaToLoCo CSVs, one <method>.csv per method")
    p.add_argument("mane_dir", type=Path,
                   help="folder with the prepare_mane.py outputs (mane_select_map_pc.tsv, "
                        "mane_select_protein.fasta, mane_cds_segments.tsv)")
    p.add_argument("--out-dir", type=Path, default=Path("."),
                   help="where to write lcr_regions.tsv, lcr_summary.txt and bed/")
    p.add_argument("--methods", nargs="+", default=DEFAULT_METHODS,
                   help="methods to process (names of the CSV files without .csv)")
    p.add_argument("--major-frac", type=float, default=0.2,
                   help="minimum fraction of a residue in an LCR to list it in major_aa")
    return p.parse_args()


def merge(intervals):
    """Merge overlapping or touching 1-based inclusive intervals."""
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def complement(merged, length):
    """Residue ranges of a protein not covered by the merged intervals."""
    out, pos = [], 1
    for s, e in merged:
        if s > pos:
            out.append((pos, s - 1))
        pos = e + 1
    if pos <= length:
        out.append((pos, length))
    return out


def describe(seq, major_frac):
    """Dominant residue, purity, residues at >= major_frac, Shannon entropy (bits)."""
    c = Counter(seq)
    n = len(seq)
    dom, k = c.most_common(1)[0]
    major = "".join(sorted(a for a, x in c.items() if x / n >= major_frac))
    ent = -sum(x / n * math.log2(x / n) for x in c.values())
    return dom, round(k / n, 3), major, abs(round(ent, 3))


def bed_lines(segs, tx, ranges, name_fn):
    """BED rows (chrom, start, end, name, score, strand) for protein ranges of transcript tx."""
    lines = []
    for i, (s, e) in enumerate(ranges):
        for chrom, bs, be, strand in protein_to_genome(segs, tx, s, e):
            lines.append((chrom, bs, be, name_fn(i, s, e), 0, strand))
    return lines


def main():
    args = parse_args()

    # Analysed MANE proteins: IDs, gene symbols, sequences and CDS segments
    m = pd.read_csv(args.mane_dir / "mane_select_map_pc.tsv", sep="\t")
    p2tx = dict(zip(m["Ensembl_prot"], m["Ensembl_nuc"]))
    p2sym = dict(zip(m["Ensembl_prot"], m["symbol"]))
    prot = {r.id: str(r.seq) for r in SeqIO.parse(args.mane_dir / "mane_select_protein.fasta", "fasta")}
    segs = load_segments(args.mane_dir / "mane_cds_segments.tsv")

    bed_dir = args.out_dir / "bed"
    bed_dir.mkdir(parents=True, exist_ok=True)
    rows, summary = [], []
    for meth in args.methods:
        t = pd.read_csv(args.results_dir / f"{meth}.csv")
        unknown = set(t["protein_id"]) - set(prot)
        if unknown:
            print(f"{meth}: {len(unknown)} protein IDs not in MANE set, skipped")
        t = t[t["protein_id"].isin(prot)].copy()

        # Coordinates must be 1-based and inside the protein. start = 0 means
        # SEG results still in segmasker's 0-based form: every SEG region would
        # be one residue off, so stop instead of silently using them.
        plen = t["protein_id"].map(lambda p: len(prot[p]))
        bad = (t["start"] < 1) | (t["end"] > plen) | (t["start"] > t["end"])
        if bad.any():
            raise SystemExit(
                f"{meth}: {bad.sum()} regions outside the protein (e.g. start = 0).\n"
                f"If these are SEG results from an older run, correct them first:\n"
                f"    python3 fix_seg_coordinates.py {args.results_dir}")

        lcr_bed, non_bed = [], []
        # protein -> merged LCR intervals
        by_prot = t.groupby("protein_id")[["start", "end"]].apply(
            lambda g: merge(zip(g["start"], g["end"]))).to_dict()

        for pid, seq in prot.items():
            merged = by_prot.get(pid, [])
            tx = p2tx[pid]
            for s, e in merged:
                sub = seq[s - 1:e]
                dom, pur, major, ent = describe(sub, args.major_frac)
                lcr_id = f"{meth}|{p2sym[pid]}|{s}-{e}"
                iv = protein_to_genome(segs, tx, s, e)
                rows.append({"lcr_id": lcr_id, "method": meth, "gene": p2sym[pid],
                             "protein": pid, "transcript": tx,
                             "start": s, "end": e, "length": e - s + 1,
                             "dominant_aa": dom, "purity": pur, "major_aa": major,
                             "entropy": ent, "chrom": iv[0][0], "strand": iv[0][3],
                             "genomic": ";".join(f"{c}:{a + 1}-{b}" for c, a, b, _ in iv),
                             "sequence": sub})
            lcr_bed += bed_lines(segs, tx, merged,
                                 lambda i, s, e, p=pid: f"{meth}|{p2sym[p]}|{s}-{e}")
            non_bed += bed_lines(segs, tx, complement(merged, len(seq)),
                                 lambda i, s, e, p=pid: f"non|{p2sym[p]}|{s}-{e}")

        for name, lines in (("lcr", lcr_bed), ("nonlcr", non_bed)):
            pd.DataFrame(lines).sort_values([0, 1]).to_csv(
                bed_dir / f"{name}_{meth}.bed", sep="\t", header=False, index=False)

        n_res = sum(e - s + 1 for v in by_prot.values() for s, e in v)
        tot = sum(len(s) for s in prot.values())
        summary.append(f"{meth:17s} {sum(len(v) for v in by_prot.values()):7d} merged LCRs "
                       f"(from {len(t)} raw) in {len(by_prot):6d} proteins, "
                       f"{n_res / tot:.2%} of residues")

    df = pd.DataFrame(rows)
    df.to_csv(args.out_dir / "lcr_regions.tsv", sep="\t", index=False)

    # Positive controls: each known repeat should be found by at least one method
    summary.append("\nPositive controls (known disease repeats):")
    for gene, aa in CONTROLS.items():
        # Repeat residue among the major residues, not necessarily the dominant
        # one: e.g. in HTT the polyQ is merged with the adjacent polyP (P > Q)
        hits = df[(df["gene"] == gene) & df["major_aa"].str.contains(aa)]
        found = ", ".join(f"{r.method} {r.start}-{r.end}" for r in hits.itertuples()) or "NOT FOUND"
        summary.append(f"  {gene:7s} poly{aa}: {found}")
    text = f"Analysed proteins: {len(prot)}\n\n" + "\n".join(summary) + "\n"
    (args.out_dir / "lcr_summary.txt").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
