#!/usr/bin/env python3
"""
Check that each MANE Select coding sequence translates exactly into its protein.

Inputs (in the current folder):
  mane_select_map_pc.tsv     protein-coding MANE Select mapping table
  mane_select_protein.fasta  Select protein sequences (Ensembl protein IDs)
  mane_select_cds.fasta      Select coding sequences (Ensembl transcript IDs)
Output:
  mane_cds_check.tsv         one row per transcript with a status
"""
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

MAP, PROT, CDS, OUT = ("mane_select_map_pc.tsv", "mane_select_protein.fasta",
                       "mane_select_cds.fasta", "mane_cds_check.tsv")

m = pd.read_csv(MAP, sep="\t")
sym_col = "symbol" if "symbol" in m.columns else None

prot = {r.id: str(r.seq).upper().rstrip("*") for r in SeqIO.parse(PROT, "fasta")}
cds = {r.id: str(r.seq).upper() for r in SeqIO.parse(CDS, "fasta")}


def first_difference(a, b):
    """1-based position of the first differing residue, or None if identical."""
    for i, (x, y) in enumerate(zip(a, b), start=1):
        if x != y:
            return i
    return None if len(a) == len(b) else min(len(a), len(b)) + 1


rows = []
for _, r in m.iterrows():
    tx, p, chrom = r["Ensembl_nuc"], r["Ensembl_prot"], r["chrom"]
    gene = r[sym_col] if sym_col else ""
    row = {"gene": gene, "transcript": tx, "protein": p, "chrom": chrom,
           "cds_len": None, "prot_len": None, "first_diff": None}

    if tx not in cds:
        rows.append({**row, "status": "cds_missing"}); continue
    if p not in prot:
        rows.append({**row, "status": "protein_missing"}); continue

    s, ref = cds[tx], prot[p]
    table = 2 if chrom == "chrM" else 1          # mitochondrial genetic code on chrM
    aa = str(Seq(s[: len(s) // 3 * 3]).translate(table=table))
    if aa.endswith("*"):                          # drop the final stop codon
        aa = aa[:-1]
    row.update(cds_len=len(s), prot_len=len(ref))

    flags = []
    if len(s) % 3:
        flags.append("length_not_multiple_of_3")

    if aa == ref:
        status = "ok"
    elif "U" in ref and aa.replace("*", "U") == ref:
        status = "ok_selenoprotein"               # UGA read as selenocysteine
    else:
        if "*" in aa:
            flags.append("internal_stop")
        flags.append("mismatch")
        row["first_diff"] = first_difference(aa, ref)
        status = ";".join(flags)

    if status.startswith("ok") and flags:
        status += ";" + ";".join(flags)
    rows.append({**row, "status": status})

res = pd.DataFrame(rows)
res.to_csv(OUT, sep="\t", index=False)

print(f"Checked {len(res)} transcripts -> {OUT}\n")
print(res["status"].value_counts().to_string())

problems = res[~res["status"].isin(["ok", "ok_selenoprotein", "cds_missing"])]
if len(problems):
    print(f"\n{len(problems)} transcripts need a look:")
    print(problems.head(30).to_string(index=False))


"""status
ok                             19167
cds_missing                       64
mismatch                          36
ok_selenoprotein                  24
internal_stop;mismatch             1
ok;length_not_multiple_of_3        1

38 transcripts need a look:
    gene         transcript           protein chrom  cds_len  prot_len  first_diff                 status
    BAG1  ENST00000634734.3 ENSP00000489189.2  chr9   1038.0     345.0         1.0               mismatch
  EIF4G2 ENST00000339995.11 ENSP00000340281.6 chr11   2724.0     907.0         1.0               mismatch
   GTF3A ENST00000381140.10 ENSP00000370532.5 chr13   1098.0     365.0         1.0               mismatch
     HCK  ENST00000375852.5 ENSP00000365012.3 chr20   1581.0     526.0         1.0               mismatch
     MYC  ENST00000621592.8 ENSP00000478887.2  chr8   1365.0     454.0         1.0               mismatch
   TEAD1  ENST00000527636.7 ENSP00000435233.2 chr11   1281.0     426.0         1.0               mismatch
   TEAD4  ENST00000359864.8 ENSP00000352926.3 chr12   1305.0     434.0         1.0               mismatch
   TEAD3 ENST00000338863.13 ENSP00000345772.8  chr6   1308.0     435.0         1.0               mismatch
   VEGFA  ENST00000672860.3 ENSP00000500082.3  chr6   1188.0     395.0         1.0               mismatch
     WT1 ENST00000452863.10 ENSP00000415516.5 chr11   1569.0     522.0         1.0               mismatch
   NR1I2  ENST00000393716.8 ENSP00000377319.3  chr3   1305.0     434.0         1.0               mismatch
   DDX17  ENST00000403230.3 ENSP00000385536.2 chr22   2190.0     729.0         1.0               mismatch
   NPTXR  ENST00000333039.4 ENSP00000327545.3 chr22   1503.0     500.0         1.0               mismatch
   NELFB  ENST00000343053.6 ENSP00000339495.6  chr9   1887.0     628.0         1.0               mismatch
SLC22A17 ENST00000354772.10 ENSP00000346824.5 chr14   1896.0     631.0         1.0               mismatch
    PIGX  ENST00000392391.9 ENSP00000376192.4  chr3    777.0     258.0         1.0               mismatch
    TUG1  ENST00000644773.3 ENSP00000501406.2 chr22    462.0     153.0         1.0               mismatch
   TRPV6  ENST00000359396.9 ENSP00000352358.5  chr7   2298.0     765.0         1.0               mismatch
   MRTFA ENST00000355630.10 ENSP00000347847.5 chr22   3096.0    1031.0         1.0               mismatch
   STIM2  ENST00000467087.7 ENSP00000419073.2  chr4   2241.0     746.0         1.0               mismatch
   HMHB1  ENST00000289448.4 ENSP00000289448.3  chr5    126.0      41.0         1.0               mismatch
  CACNG8  ENST00000270458.4 ENSP00000270458.3 chr19   1278.0     425.0         1.0               mismatch
   PANK2  ENST00000610179.7 ENSP00000477429.2 chr20   1383.0     460.0         1.0               mismatch
   WDR26  ENST00000414423.9 ENSP00000408108.4  chr1   2286.0     761.0         1.0               mismatch
  NFKBID  ENST00000641389.3 ENSP00000493265.2 chr19   1368.0     455.0         1.0               mismatch
  TXNRD3  ENST00000524230.9 ENSP00000430031.4  chr3   1932.0     643.0         1.0 internal_stop;mismatch
  KCTD11  ENST00000333751.8 ENSP00000328352.5 chr17    816.0     271.0         1.0               mismatch
  RNF187  ENST00000305943.9 ENSP00000306396.9  chr1    708.0     235.0         1.0               mismatch
C16orf82  ENST00000505035.3 ENSP00000500329.2 chr16    465.0     154.0         1.0               mismatch
  R3HCC1 ENST00000265806.12 ENSP00000265806.8  chr8   1323.0     440.0         1.0               mismatch"""