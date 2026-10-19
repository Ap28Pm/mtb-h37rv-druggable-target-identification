# filter_short_proteins.py
input_faa = "output_faa/core_high_confidence.faa"
output_faa = "output_faa/core_filtered.faa"

min_length = 100
kept, removed = 0, 0

with open(input_faa) as f, open(output_faa, "w") as out:
    seq_id, seq = "", ""
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            # save previous sequence if valid
            if seq_id:
                if len(seq) >= min_length:
                    out.write(f">{seq_id}\n{seq}\n")
                    kept += 1
                else:
                    removed += 1
            seq_id = line[1:].split()[0]
            seq = ""
        else:
            seq += line
    # check last sequence
    if seq_id:
        if len(seq) >= min_length:
            out.write(f">{seq_id}\n{seq}\n")
            kept += 1
        else:
            removed += 1

print(f"Proteins kept: {kept}")
print(f"Proteins removed (<{min_length} aa): {removed}")
print(f"Filtered file written to: {output_faa}")
