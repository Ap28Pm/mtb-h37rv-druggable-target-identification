#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import csv

# Input / Output files
input_faa = "output_faa/CNHEVCME.faa"
output_csv = "CNHEVCME_characterization.csv"
filtered_faa = "output_faa/CNHEVCME_characterized.faa"

# Filter criteria
min_mw = 10000       # 10 kDa
max_mw = 110000      # 110 kDa
max_instability = 40.0
max_gravy = 0.5      # optional

# Prepare CSV
with open(output_csv, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow([
        "Accession", "Sequence_Length", "Molecular_Weight(Da)", 
        "Instability_Index", "GRAVY", "Aliphatic_Index", "Pass_Filter"
    ])

    kept_records = []

    # Parse FASTA
    for record in SeqIO.parse(input_faa, "fasta"):
        seq = str(record.seq)
        prot = ProteinAnalysis(seq)
        length = len(seq)  # <- defined before using in calculations
        mw = prot.molecular_weight()
        instability = prot.instability_index()
        gravy = prot.gravy()

        # Calculate aliphatic index manually
        aa_counts = prot.count_amino_acids()
        aliphatic = (
            (aa_counts.get("A", 0) / length * 100) +
            2.9 * (aa_counts.get("V", 0) / length * 100) +
            3.9 * ((aa_counts.get("I", 0) / length * 100) + (aa_counts.get("L", 0) / length * 100))
        )

        # Apply filtering criteria
        pass_filter = (min_mw <= mw <= max_mw) and (instability < max_instability) and (gravy <= max_gravy)

        # Write characterization to CSV
        writer.writerow([
            record.id, length, round(mw, 2), round(instability, 2), 
            round(gravy, 3), round(aliphatic, 2), pass_filter
        ])

        # Keep passing sequences for new FASTA
        if pass_filter:
            kept_records.append(record)

# Write filtered FASTA
with open(filtered_faa, "w") as out_f:
    SeqIO.write(kept_records, out_f, "fasta")

# Summary
total_sequences = len(list(SeqIO.parse(input_faa, "fasta")))
print(f"Total sequences analyzed: {total_sequences}")
print(f"Sequences passing filter: {len(kept_records)}")
print(f"Characterization CSV saved as: {output_csv}")
print(f"Filtered FASTA saved as: {filtered_faa}")
