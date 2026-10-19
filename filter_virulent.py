#!/usr/bin/env python3
import csv
from Bio import SeqIO

# Input files
blast_file = "blast_results/CNHE_vs_V.tsv"  # BLASTp results
query_fasta = "output_faa/CNHE.faa"         # Essential non-human proteins
output_fasta = "output_faa/CNHEV.faa"

# Parse BLAST results to get all query sequences with any hit
virulent_queries = set()
with open(blast_file) as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        if len(row) < 1:
            continue
        qseqid = row[0]
        virulent_queries.add(qseqid)

# Write sequences that have hits in the NHEV DB
with open(output_fasta, "w") as out_f:
    for record in SeqIO.parse(query_fasta, "fasta"):
        if record.id in virulent_queries:
            SeqIO.write(record, out_f, "fasta")

print(f"Total virulent proteins identified: {len(virulent_queries)}")
print(f"FASTA written to: {output_fasta}")
