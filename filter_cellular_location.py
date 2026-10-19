import csv
from Bio import SeqIO
from collections import defaultdict

def filter_psortb_localizations(psortb_file: str, fasta_file: str, 
                                output_fasta: str, per_loc: bool = True):
    """
    Filters a FASTA file based on PSORTb final localization results.

    Args:
        psortb_file (str): Path to PSORTb results (.txt) file.
        fasta_file (str): Path to the protein FASTA file.
        output_fasta (str): Path to write combined filtered sequences.
        per_loc (bool): If True, also generate separate FASTAs per localization.
    """
    target_localizations = {"cellwall", "cytoplasmicmembrane", "extracellular"}

    # Read PSORTb results
    localization_map = defaultdict(list)  # localization -> list of seq IDs
    total_proteins = 0
    with open(psortb_file, 'r', newline='') as pf:
        reader = csv.DictReader(pf, delimiter='\t')
        if 'SeqID' not in reader.fieldnames or 'Final_Localization' not in reader.fieldnames:
            raise ValueError("PSORTb results file must have 'SeqID' and 'Final_Localization' columns.")
        
        for row in reader:
            seq_id_full = row['SeqID'].split()[0]  # take first token
            final_loc = row['Final_Localization'].strip().lower()
            total_proteins += 1

            # Split multi-localization (PSORTb separates by ';' or ',')
            locs = [loc.strip() for loc in final_loc.replace(',', ';').split(';') if loc.strip()]
            for loc in locs:
                if loc in target_localizations:
                    localization_map[loc].append(seq_id_full)

    # Print summary
    print("### PSORTb Localization Summary ###")
    print(f"Total proteins analyzed: {total_proteins}\n")
    print(f"{'Localization':25s} {'Count':>6s} {'Percentage':>10s}")
    for loc in target_localizations:
        count = len(localization_map[loc])
        percentage = (count / total_proteins) * 100
        print(f"{loc:25s} {count:6d} {percentage:10.2f}%")
    
    # Load sequences and filter
    seq_records = list(SeqIO.parse(fasta_file, "fasta"))
    combined_filtered = []
    per_loc_records = defaultdict(list)
    
    seq_id_set = set(seq.id for seq in seq_records)
    for loc, ids in localization_map.items():
        for seq_id in ids:
            if seq_id in seq_id_set:
                record = next(rec for rec in seq_records if rec.id == seq_id)
                combined_filtered.append(record)
                per_loc_records[loc].append(record)

    # Remove duplicates in combined list
    combined_filtered = list({rec.id: rec for rec in combined_filtered}.values())

    # Write combined FASTA
    SeqIO.write(combined_filtered, output_fasta, "fasta")
    print(f"\nCombined filtered sequences written to: {output_fasta}")

    # Write separate FASTAs if requested
    if per_loc:
        for loc, records in per_loc_records.items():
            outfile = f"{loc}_proteins.faa"
            SeqIO.write(records, outfile, "fasta")
            print(f"{len(records)} sequences written to {outfile}")

# ------------------ Usage ------------------
if __name__ == "__main__":
    psortb_file = "PSORTb_result.txt"           # Replace with your file
    input_fasta = "output_faa/CNHEV.faa"  # Replace with your FASTA
    output_fasta = "output_faa/CHNEV_cellular_location.faa"

    filter_psortb_localizations(psortb_file, input_fasta, output_fasta)
