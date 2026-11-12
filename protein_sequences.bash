col=$(head -1 mane_select_map_pc.tsv | tr '\t' '\n' | grep -n -x "Ensembl_prot" | cut -d: -f1)
cut -f$col mane_select_map_pc.tsv | tail -n +2 > select_ensp.txt
seqkit grep -f select_ensp.txt MANE.GRCh38.v1.5.ensembl_protein.faa.gz > mane_select_protein.fasta
seqkit stats mane_select_protein.fasta