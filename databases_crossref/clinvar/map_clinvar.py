"""
Map ClinVar variants onto MANE Select proteins and flag those inside LCRs.

Usage:
    python3 map_clinvar.py ~/data/clinvar ~/data/mane lcr_regions.tsv

Inputs:
  clinvar/variant_summary.txt.gz   classification, review status, origin
  clinvar/hgvs4variation.txt.gz    protein-level HGVS of each variant
  mane/mane_select_map_pc.tsv      RefSeq protein -> Ensembl protein of MANE Select
  mane/mane_select_protein.fasta   protein sequences (reference residue check)
  lcr_regions.tsv                  merged LCRs per method (build_lcr_table.py)

What it does:
  - keeps protein-level HGVS on the MANE Select RefSeq protein of each gene
    (version-independent), so every variant has one position on the protein
    used for LCR detection
  - parses the change: missense, nonsense, synonymous, frameshift, in-frame
    deletion/insertion/duplication/delins, repeat-length change, start loss,
    stop loss; variants without a protein change (intronic, UTR) are not
    included
  - checks that the reference residue in the HGVS matches the MANE sequence
    (ref_ok); mismatches come from older protein versions and should be dropped
  - adds classification group (P/LP, B/LB, VUS, conflicting, other), review
    stars (0-4) and origin
  - for each LCR method, flags whether the variant overlaps an LCR and which one

Outputs:
  clinvar_mane_variants.tsv  one row per variant per MANE protein
  clinvar_summary.txt        counts, and share of P/LP variants inside LCRs
                             compared with the share of residues inside LCRs
"""
import gzip
import re
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import seq1

clinvar, mane, lcr_file = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
METHODS = ["SEG_strict", "SEG_intermediate", "CAST", "fLPS_strict"]


def read_clinvar_table(path, usecols_wanted, chunk_filter=None):
    """Read a ClinVar tab-delimited file whose header line starts with '#'.
    Comment lines before the header are skipped."""
    with gzip.open(path, "rt") as fh:
        skip = 0
        for line in fh:
            if line.startswith("#") and not line.startswith("##") and "\t" in line:
                break
            skip += 1
    header = pd.read_csv(path, sep="\t", skiprows=skip, nrows=0).columns
    names = [c.lstrip("#") for c in header]
    cols = [c for c in usecols_wanted if c in names]
    parts = []
    for chunk in pd.read_csv(path, sep="\t", skiprows=skip + 1, header=None,
                             names=names, usecols=cols, dtype=str,
                             chunksize=500_000):
        parts.append(chunk_filter(chunk) if chunk_filter else chunk)
    return pd.concat(parts, ignore_index=True), names


def pick(names, candidates, what):
    for c in candidates:
        if c in names:
            return c
    sys.exit(f"No column for {what} in ClinVar file. Columns are:\n{names}")


# --- MANE --------------------------------------------------------------------
m = pd.read_csv(mane / "mane_select_map_pc.tsv", sep="\t")
m["np"] = m["RefSeq_prot"].str.split(".").str[0]
np2ensp = dict(zip(m["np"], m["Ensembl_prot"]))
ensp2sym = dict(zip(m["Ensembl_prot"], m["symbol"]))
prot = {r.id: str(r.seq) for r in SeqIO.parse(mane / "mane_select_protein.fasta", "fasta")}

# --- ClinVar protein HGVS ----------------------------------------------------
H_COLS = ["VariationID", "ProteinExpression"]
hg, hnames = read_clinvar_table(
    clinvar / "hgvs4variation.txt.gz", H_COLS,
    lambda c: c[c["ProteinExpression"].str.startswith("NP_", na=False)])
hg["np"] = hg["ProteinExpression"].str.split(".").str[0]
hg = hg[hg["np"].isin(np2ensp)].drop_duplicates(["VariationID", "ProteinExpression"])
hg["protein"] = hg["np"].map(np2ensp)
hg["hgvs_p"] = hg["ProteinExpression"].str.split(":").str[1]
print(f"Protein HGVS on MANE Select proteins: {len(hg)}")

# --- ClinVar classifications -------------------------------------------------
with gzip.open(clinvar / "variant_summary.txt.gz", "rt") as fh:
    vnames = fh.readline().rstrip("\n").lstrip("#").split("\t")
col_sig = pick(vnames, ["ClinicalSignificance", "GermlineClassification"], "classification")
col_rev = pick(vnames, ["ReviewStatus", "GermlineReviewStatus"], "review status")
col_ori = pick(vnames, ["OriginSimple", "Origin"], "origin")
vs, _ = read_clinvar_table(
    clinvar / "variant_summary.txt.gz",
    ["VariationID", "Type", "Assembly", "PhenotypeList", col_sig, col_rev, col_ori],
    lambda c: c[c["Assembly"] == "GRCh38"])
vs = vs.drop_duplicates("VariationID").rename(
    columns={col_sig: "clinical_significance", col_rev: "review_status",
             col_ori: "origin", "Type": "variant_type", "PhenotypeList": "conditions"})


def group(sig):
    s = str(sig).lower()
    if "conflicting" in s:
        return "conflicting"
    if "pathogenic" in s and "benign" not in s:
        return "P/LP"
    if "benign" in s and "pathogenic" not in s:
        return "B/LB"
    if "uncertain" in s:
        return "VUS"
    return "other"


def stars(rev):
    r = str(rev).lower()
    if "practice guideline" in r:
        return 4
    if "expert panel" in r:
        return 3
    if "multiple submitters" in r and "no conflict" in r:
        return 2
    if "criteria provided" in r:
        return 1
    return 0


vs["class_group"] = vs["clinical_significance"].map(group)
vs["stars"] = vs["review_status"].map(stars)
v = hg.merge(vs.drop(columns="Assembly"), on="VariationID", how="inner")

# --- Parse protein HGVS ------------------------------------------------------
AA = r"(?:[A-Z][a-z]{2}|\*)"
RX = re.compile(rf"^({AA})(\d+)(?:_({AA})(\d+))?(.*)$")


def parse(hgvs):
    """-> (consequence, aa_start, aa_end, ref_first, ref_last)"""
    p = hgvs[2:] if hgvs.startswith("p.") else hgvs
    p = p.strip("()")
    if p in ("?", "=", "0", ""):
        return ("unknown" if p != "=" else "synonymous", None, None, None, None)
    mt = RX.match(p)
    if not mt:
        return ("unparsed", None, None, None, None)
    r1, s, r2, e, rest = mt.groups()
    s = int(s)
    e = int(e) if e else s
    if rest.startswith("["):
        cons = "repeat_change"
    elif "fs" in rest:
        cons = "frameshift"
    elif rest.startswith("delins"):
        cons = "delins"
    elif rest.startswith("del"):
        cons = "inframe_deletion"
    elif rest.startswith("dup"):
        cons = "duplication"
    elif rest.startswith("ins"):
        cons = "inframe_insertion"
    elif "ext" in rest or r1 in ("Ter", "*"):
        cons = "stop_loss"
    elif rest == "=":
        cons = "synonymous"
    elif rest in ("Ter", "*"):
        cons = "nonsense"
    elif s == 1 and r1 == "Met":
        cons = "start_loss"
    elif re.fullmatch(r"[A-Z][a-z]{2}", rest):
        cons = "missense"
    else:
        cons = "other"
    return (cons, s, e, r1, r2 or r1)


parsed = v["hgvs_p"].map(parse)
v[["consequence", "aa_start", "aa_end", "ref_first", "ref_last"]] = pd.DataFrame(
    parsed.tolist(), index=v.index)
v = v[v["aa_start"].notna()].copy()
v["aa_start"] = v["aa_start"].astype(int)
v["aa_end"] = v["aa_end"].astype(int)


def ref_ok(row):
    seq = prot.get(row["protein"], "")
    for pos, res in ((row["aa_start"], row["ref_first"]), (row["aa_end"], row["ref_last"])):
        if res in ("Ter", "*"):
            if pos != len(seq) + 1:
                return False
        elif pos > len(seq) or seq[pos - 1] != seq1(res):
            return False
    return True


v["ref_ok"] = v.apply(ref_ok, axis=1)
v["gene"] = v["protein"].map(ensp2sym)

# --- Overlap with LCRs -------------------------------------------------------
lcr = pd.read_csv(lcr_file, sep="\t")
for meth in METHODS:
    sub = lcr[lcr["method"] == meth]
    by_prot = {p: list(zip(g["start"], g["end"], g["lcr_id"])) for p, g in sub.groupby("protein")}

    def hit(row, by_prot=by_prot):
        for s, e, lid in by_prot.get(row["protein"], []):
            if row["aa_start"] <= e and row["aa_end"] >= s:
                return lid
        return ""

    v[f"lcr_{meth}"] = v.apply(hit, axis=1)

cols = ["VariationID", "gene", "protein", "ProteinExpression", "hgvs_p", "consequence",
        "aa_start", "aa_end", "ref_ok", "variant_type", "clinical_significance",
        "class_group", "review_status", "stars", "origin", "conditions"] \
       + [f"lcr_{x}" for x in METHODS]
v[cols].to_csv("clinvar_mane_variants.tsv", sep="\t", index=False)

# --- Summary -----------------------------------------------------------------
ok = v[v["ref_ok"]]
lines = [f"Variants with protein position on MANE Select: {len(v)}",
         f"  reference residue matches MANE: {len(ok)} "
         f"({len(v) - len(ok)} mismatches - drop these)",
         "", "Consequences (ref_ok only):", ok["consequence"].value_counts().to_string(),
         "", "Classification groups (ref_ok only):", ok["class_group"].value_counts().to_string(),
         "", "First look - P/LP missense, >=1 star, ref_ok:"]
plp = ok[(ok["class_group"] == "P/LP") & (ok["stars"] >= 1) & (ok["consequence"] == "missense")]
blb = ok[(ok["class_group"] == "B/LB") & (ok["stars"] >= 1) & (ok["consequence"] == "missense")]
total_res = sum(len(s) for s in prot.values())
for meth in METHODS:
    sub = lcr[lcr["method"] == meth]
    res_share = sub["length"].sum() / total_res
    p_in = (plp[f"lcr_{meth}"] != "").mean() if len(plp) else 0
    b_in = (blb[f"lcr_{meth}"] != "").mean() if len(blb) else 0
    lines.append(f"  {meth:17s} residues in LCRs {res_share:6.2%} | "
                 f"P/LP in LCRs {p_in:6.2%} (n={len(plp)}) | B/LB in LCRs {b_in:6.2%} (n={len(blb)})")
lines += ["", "This is a crude first look across all proteins: the real comparison is",
          "within proteins (LCR vs rest of the same protein) with a mutability model."]
text = "\n".join(lines) + "\n"
Path("clinvar_summary.txt").write_text(text)
print(text)