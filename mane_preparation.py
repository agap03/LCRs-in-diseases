import pandas as pd

s = pd.read_csv("MANE.GRCh38.v1.5.summary.txt.gz", sep="\t")
s.columns = s.columns.str.lstrip("#")
print(s.columns.tolist())
print(s["MANE_status"].value_counts())

sel = s[s["MANE_status"] == "MANE Select"].copy()

# GRCh38_chr holds RefSeq accessions like NC_000001.11 -> convert to chr1
def to_chr(acc):
    n = int(acc.split("_")[1].split(".")[0])
    return {23: "chrX", 24: "chrY", 12920: "chrM"}.get(n, f"chr{n}")
sel["chrom"] = sel["GRCh38_chr"].map(to_chr)

pc = sel[sel["Ensembl_prot"].notna()]
pc.to_csv("mane_select_map_pc.tsv", sep="\t", index=False)