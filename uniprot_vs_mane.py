from Bio import SeqIO
import gzip
import pandas as pd

uni = {}
for r in SeqIO.parse(gzip.open("../uniprot/UP000005640_9606.fasta.gz", "rt"), "fasta"):
    if r.id.startswith("sp|"):
        uni[str(r.seq)] = r.id.split("|")[1]

rows = []
for r in SeqIO.parse("mane_select_protein.fasta", "fasta"):
    rows.append((r.id, uni.get(str(r.seq).rstrip("*"))))

m = pd.DataFrame(rows, columns=["ensp", "uniprot_identical"])
m.to_csv("mane_uniprot_identity.tsv", sep="\t", index=False)
print(m["uniprot_identical"].notna().mean())

# result:
# 0.9109003265433059 - the fraction of genes where UniProt and MANE positions are interchangeable. 