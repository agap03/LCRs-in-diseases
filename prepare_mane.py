"""
Prepare MANE Select v1.5 for the analysis. Run from ~/data/mane.

Inputs (in the current directory unless noted):
  MANE.GRCh38.v1.5.summary.txt.gz          one row per MANE transcript (IDs, status, chromosome)
  MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz  genomic annotation (genes, exons, CDS, ...)
  MANE.GRCh38.v1.5.ensembl_protein.faa.gz  protein sequences (Ensembl ENSP IDs)
  ../uniprot/UP000005640_9606.fasta.gz     UniProt human reference proteome

Outputs:
  mane_select_map_pc.tsv     analysed genes: protein-coding Select genes on primary
                             chromosomes, with RefSeq/Ensembl IDs and chrom
  mane_excluded_genes.txt    excluded transcripts with reasons
  mane_primary.gtf           MANE annotation on chr1-22, X, Y, M
  mane_select_protein.fasta  protein sequences of the analysed genes
  mane_cds_segments.tsv      CDS segments of each transcript (GTF coordinates)
  mane_uniprot_identity.tsv  identical Swiss-Prot entry for each protein, if any
"""
import gzip
import re
from datetime import date
import pandas as pd
from Bio import SeqIO

PREFIX = "MANE.GRCh38.v1.5"                          # common start of all MANE file names
UNIPROT = "../uniprot/UP000005640_9606.fasta.gz"     # human reference proteome
# Primary chromosomes: chr1-22, chrX, chrY, chrM. Anything else (e.g. fix-patch
# contigs such as "chr1_KN538360v1_fix") is not part of the hg38 primary assembly.
PRIMARY = re.compile(r"^chr([0-9]+|X|Y|M)$")


def to_chr(acc):
    """NC_000001.11 -> chr1; other accessions (patches) unchanged.

    The MANE summary gives chromosomes as RefSeq accessions. NC_0000NN are the
    nuclear chromosomes (23 = X, 24 = Y) and NC_012920 is the mitochondrion.
    """
    if acc.startswith("NC_012920"):
        return "chrM"
    if not acc.startswith("NC_0000"):
        return acc
    n = int(acc[3:].split(".")[0])          # "NC_000001.11" -> "000001" -> 1
    return {23: "chrX", 24: "chrY"}.get(n, f"chr{n}")


def main():
    # 1. MANE Select genes from the summary table
    s = pd.read_csv(f"{PREFIX}.summary.txt.gz", sep="\t")
    s.columns = s.columns.str.lstrip("#")       # first header starts with "#"
    # Ensembl transcript ID -> "MANE Select" / "MANE Plus Clinical" (used in exclusions)
    status = dict(zip(s["Ensembl_nuc"], s["MANE_status"]))
    # Only MANE Select (one transcript per gene); Plus Clinical transcripts are extra
    sel = s[s["MANE_status"] == "MANE Select"].copy()
    sel["chrom"] = sel["GRCh38_chr"].map(to_chr)

    # 2. Keep primary chromosomes; collect CDS segments and transcript info
    #    The GTF is read line by line: lines on primary chromosomes are copied to
    #    mane_primary.gtf, and every CDS line is stored as a segment.
    #    GTF columns: 0 seqname, 2 feature, 3 start, 4 end, 6 strand, 8 attributes.
    tx_re, gene_re = re.compile(r'transcript_id "([^"]+)"'), re.compile(r'gene_name "([^"]+)"')
    tx_info = {}    # transcript -> (gene name, contig), for every transcript in the GTF
    cds_rows = []   # (transcript, chrom, start, end, strand) for each CDS line
    with gzip.open(f"{PREFIX}.ensembl_genomic.gtf.gz", "rt") as fin, open("mane_primary.gtf", "w") as fout:
        for line in fin:
            if line.startswith("#"):       # header/comment lines: keep as they are
                fout.write(line)
                continue
            f = line.rstrip("\n").split("\t")

            m = tx_re.search(f[8])
            if m is None:                  # gene-level lines have no transcript_id
                if PRIMARY.match(f[0]):
                    fout.write(line)       # keep them in mane_primary.gtf
                continue
            tx = m.group(1)

            # Remember gene name and contig the first time a transcript is seen
            if tx not in tx_info:
                g = gene_re.search(f[8])
                tx_info[tx] = (g.group(1) if g else "", f[0])
            if PRIMARY.match(f[0]):
                fout.write(line)
                if f[2] == "CDS":
                    cds_rows.append((tx, f[0], int(f[3]), int(f[4]), f[6]))
    # Transcripts with at least one CDS on a primary chromosome = protein-coding
    has_cds = {r[0] for r in cds_rows}

    # 3. Exclusions: fix-patch transcripts and non-coding transcripts
    #    Written to a file so it is documented which transcripts were dropped and why.
    excl = []
    for tx, (gene, contig) in tx_info.items():
        if not PRIMARY.match(contig):
            excl.append(("fix_patch_contig", gene, tx, status.get(tx, ""), contig))
        elif tx not in has_cds:
            excl.append(("non_coding", gene, tx, status.get(tx, ""), contig))
    excl = pd.DataFrame(excl, columns=["reason", "gene", "transcript", "mane_status", "contig"])
    excl = excl.sort_values(["reason", "gene"])
    counts = excl.groupby(["reason", "mane_status"]).size()   # numbers for the file header
    with open("mane_excluded_genes.txt", "w") as out:
        # Header: date, meaning of each reason, counts per reason and MANE status
        out.write(f"# Transcripts excluded from {PREFIX} ({date.today()})\n"
                "# fix_patch_contig: only on GRCh38 fix-patch contigs, absent from hg38 and variant databases\n"
                "# non_coding: no CDS (non-coding RNA), outside protein-level analysis\n")
        for (reason, st), n in counts.items():
            out.write(f"# {reason}\t{st}\t{n}\n")
        excl.to_csv(out, sep="\t", index=False)               # then the full table

    # 4. Analysed genes: protein-coding Select transcripts with CDS on primary chromosomes
    #    (must have an Ensembl protein ID and CDS lines kept in step 2)
    pc = sel[sel["Ensembl_prot"].notna() & sel["Ensembl_nuc"].isin(has_cds)]
    pc.to_csv("mane_select_map_pc.tsv", sep="\t", index=False)
    keep_tx, keep_prot = set(pc["Ensembl_nuc"]), set(pc["Ensembl_prot"])

    # CDS segments of the analysed transcripts only, in genomic order
    # (cds_mapping.load_segments reverses them for minus-strand genes)
    seg = pd.DataFrame(cds_rows, columns=["transcript_id", "chrom", "start", "end", "strand"])
    seg = seg[seg["transcript_id"].isin(keep_tx)].sort_values(["transcript_id", "start"])
    seg.to_csv("mane_cds_segments.tsv", sep="\t", index=False)

    # 5. Protein sequences of the analysed genes
    #    Upper-cased, with a trailing stop "*" removed if present.
    prot = {}
    with gzip.open(f"{PREFIX}.ensembl_protein.faa.gz", "rt") as fh:
        for r in SeqIO.parse(fh, "fasta"):
            if r.id in keep_prot:
                prot[r.id] = str(r.seq).upper().rstrip("*")
    with open("mane_select_protein.fasta", "w") as out:
        for p, seq in prot.items():
            out.write(f">{p}\n{seq}\n")

    # 6. Identity with reviewed UniProt sequences
    #    Build sequence -> Swiss-Prot accession (reviewed "sp|" entries only, not TrEMBL "tr|").
    #    If several entries share a sequence, the first one is kept (setdefault).
    uni = {}
    with gzip.open(UNIPROT, "rt") as fh:
        for r in SeqIO.parse(fh, "fasta"):
            if r.id.startswith("sp|"):
                uni.setdefault(str(r.seq).upper(), r.id.split("|")[1])   # "sp|P04637|P53_HUMAN" -> "P04637"
    # For each MANE protein: accession of the Swiss-Prot entry with exactly the
    # same sequence, or empty if there is none (only 100% identical matches count)
    ident = pd.DataFrame({"symbol": pc["symbol"], "Ensembl_prot": pc["Ensembl_prot"]})
    ident["uniprot_identical"] = ident["Ensembl_prot"].map(lambda p: uni.get(prot.get(p)))
    ident.to_csv("mane_uniprot_identity.tsv", sep="\t", index=False)

    # Summary on screen
    print(f"MANE Select genes:                 {len(sel)}")
    print(f"Excluded transcripts:\n{counts.to_string()}")
    print(f"Analysed protein-coding genes:     {len(pc)}")
    print(f"Proteins written:                  {len(prot)}")
    print(f"Identical to Swiss-Prot:           {ident['uniprot_identical'].notna().mean():.1%}")


if __name__ == "__main__":
    main()
