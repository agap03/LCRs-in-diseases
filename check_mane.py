"""
Build and verify MANE Select coding sequences, and test protein-to-genome
conversion. Run from ~/data/mane after prepare_mane.py.

What it checks:
  1. Coding sequences rebuilt from hg38 + the CDS segments translate to exactly
     the MANE protein sequences (so the segments and the genome agree).
  2. cds_mapping.protein_to_genome() gives the right genomic intervals: random
     protein regions are mapped to the genome, translated back and compared.

Inputs:
  ../genome/hg38.fa           reference genome (with .fai index for pyfaidx)
  mane_select_map_pc.tsv      from prepare_mane.py
  mane_select_protein.fasta   from prepare_mane.py
  mane_cds_segments.tsv       from prepare_mane.py (read by cds_mapping.load_segments)

Outputs:
  mane_select_cds.fasta  coding sequences of the analysed genes, from hg38
  mane_cds_check.tsv     translation check per transcript (gene, lengths, status)
  mane_qc_summary.txt    status counts and mapping-test result
"""
import os
import random
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from pyfaidx import Fasta

def main():
    # Make cds_mapping.py importable when the script is run from another
    # directory (~/data/mane): add the script's own folder to the import path.
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from cds_mapping import load_segments, protein_to_genome, coding_dna  # noqa: E402

    GENOME = "../genome/hg38.fa"

    genome = Fasta(GENOME)                              # indexed, random access to hg38
    segs = load_segments()                              # transcript -> CDS segments in coding order
    m = pd.read_csv("mane_select_map_pc.tsv", sep="\t") # analysed genes
    prot = {r.id: str(r.seq) for r in SeqIO.parse("mane_select_protein.fasta", "fasta")}

    # 1. Coding sequences from the genome, via the CDS segments
    #    Segments are GTF coordinates (1-based, inclusive); s - 1 converts them
    #    to the 0-based half-open intervals that coding_dna() expects.
    cds = {}
    with open("mane_select_cds.fasta", "w") as out:
        for tx in m["Ensembl_nuc"]:
            cds[tx] = coding_dna(genome, [(c, s - 1, e, st) for c, s, e, st in segs[tx]])
            out.write(f">{tx}\n{cds[tx]}\n")

    # 2. Translation check, labelling known exceptions
    #    Translate each CDS and compare with the MANE protein. Known biological
    #    cases where a plain translation differs are corrected and labelled, so
    #    only unexplained differences end up as MISMATCH.
    rows = []
    for _, r in m.iterrows():
        tx, ref = r["Ensembl_nuc"], prot[r["Ensembl_prot"]]
        s = cds[tx]
        # Translate only complete codons. Mitochondrial genes use the vertebrate
        # mitochondrial genetic code (NCBI table 2), all others the standard one (table 1).
        aa = str(Seq(s[: len(s) // 3 * 3]).translate(table=2 if r["chrom"] == "chrM" else 1))
        aa = aa[:-1] if aa.endswith("*") else aa         # drop final stop, if included
        labels = []
        if ref.startswith("M") and aa and aa[0] != "M":     # non-AUG start codon read as Met
            # e.g. a CUG start codon translates to L, but the ribosome inserts Met
            aa = "M" + aa[1:]
            labels.append("nonAUG_start")
        if "U" in ref:                                    # UGA read as selenocysteine
            # In selenoproteins some UGA codons encode selenocysteine (U) instead
            # of stop (*): put U wherever the reference has U and we have *.
            aa = "".join("U" if a == "*" and b == "U" else a for a, b in zip(aa, ref)) + aa[len(ref):]
            labels.append("selenoprotein")
        if len(s) % 3:
            # CDS length not divisible by 3 (e.g. mitochondrial genes whose stop
            # codon is completed only by polyadenylation)
            labels.append("length_not_multiple_of_3")
        status = "ok" if aa == ref else "MISMATCH"
        # status column: "ok" or "MISMATCH", followed by any labels, separated by ";"
        rows.append((r["symbol"], tx, r["chrom"], len(s), len(ref), ";".join([status] + labels)))

    res = pd.DataFrame(rows, columns=["gene", "transcript", "chrom", "cds_len", "prot_len", "status"])
    res.to_csv("mane_cds_check.tsv", sep="\t", index=False)

    # 3. Mapping test: random protein regions -> genome -> translate -> compare
    #    Uses only transcripts with status exactly "ok" (no special cases), so a
    #    plain translation must reproduce the protein region exactly.
    #    chrM is left out: its genes need the mitochondrial code (table 2), while
    #    this test translates with the standard code.
    random.seed(1)                                      # fixed seed: reproducible test
    plain = list(res.loc[(res["status"] == "ok") & (res["chrom"] != "chrM"), "transcript"])
    tx2prot = dict(zip(m["Ensembl_nuc"], m["Ensembl_prot"]))
    n_test, fails = 500, 0
    for tx in random.sample(plain, min(n_test, len(plain))):   # up to 500 different transcripts
        p = prot[tx2prot[tx]]
        # Random region a..b (1-based, inclusive), 1 to 61 residues long,
        # cut at the protein end
        a = random.randint(1, len(p))
        b = min(len(p), a + random.randint(0, 60))
        # Protein region -> genomic intervals -> DNA from hg38
        dna = coding_dna(genome, protein_to_genome(segs, tx, a, b))
        # Translate and compare with the same region of the protein
        # (True counts as 1, so fails is the number of wrong regions)
        fails += str(Seq(dna).translate()) != p[a - 1:b]

    # Write and print: count of each status, then the mapping-test result
    summary = (res["status"].value_counts().to_string()
            + f"\n\nMapping test: {min(n_test, len(plain)) - fails}/{min(n_test, len(plain))} correct\n")
    open("mane_qc_summary.txt", "w").write(summary)
    print(summary)
    # List any transcripts whose translation still differs after the corrections above
    if res["status"].str.startswith("MISMATCH").any():
        print("Unexplained mismatches:\n", res[res["status"].str.startswith("MISMATCH")].to_string(index=False))

if __name__ == "__main__":
    main()
