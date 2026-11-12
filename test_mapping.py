import random
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from pyfaidx import Fasta
from cds_mapping import load_segments, protein_to_genome

genome = Fasta("../genome/hg38.fa")
segs = load_segments()
check = pd.read_csv("mane_cds_check.tsv", sep="\t")
ok = check[check["status"] == "ok"]
tx2prot = dict(zip(ok["transcript"], ok["protein"]))
prot = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse("mane_select_protein.fasta", "fasta")}


def coding_dna(intervals):
    """Coding-strand DNA for intervals given in coding order."""
    pieces = []
    for chrom, start, end, strand in intervals:
        seq = str(genome[chrom][start:end]).upper()
        pieces.append(str(Seq(seq).reverse_complement()) if strand == "-" else seq)
    return "".join(pieces)


random.seed(1)
fails = 0
for tx in random.sample(sorted(tx2prot), 500):
    p = prot[tx2prot[tx]]
    a = random.randint(1, len(p))
    b = min(len(p), a + random.randint(0, 60))
    ivs = protein_to_genome(segs, tx, a, b)
    if str(Seq(coding_dna(ivs)).translate()) != p[a - 1:b]:
        fails += 1
        print("FAIL", tx, a, b, ivs)
print(f"{500 - fails}/500 correct")