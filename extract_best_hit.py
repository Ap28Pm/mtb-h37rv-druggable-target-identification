from Bio import SeqIO
import pandas as pd

# Input files
best_hits_csv = "blast_results/query_vs_core_best_hits.csv"
query_fasta = "H3Rv7_proteome.faa"
output_faa = "output_faa/core_high_confidence_429.faa"

# Load best-hit IDs (high_match == True)
df = pd.read_csv(best_hits_csv)
high_match_ids = set(df[df['high_match'] == True]['qseqid'].tolist())

print(f"Number of high-confidence matches: {len(high_match_ids)}")

# Parse query fasta and extract only high-confidence sequences
sequences = []
for record in SeqIO.parse(query_fasta, "fasta"):
    if record.id in high_match_ids:
        sequences.append(record)

# Write to output FASTA
SeqIO.write(sequences, output_faa, "fasta")
print(f"{len(sequences)} high-confidence core sequences written to {output_faa}")

