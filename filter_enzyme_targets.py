import pandas as pd
from Bio import SeqIO
import os

def extract_enzyme_targets(interpro_file, input_fasta, output_fasta):
    """
    Identifies enzyme proteins from InterProScan results using the corrected column index (12), 
    extracts their sequences from the input FASTA, and saves them to a new FASTA file.

    Args:
        interpro_file (str): Path to the InterProScan TSV results file.
        input_fasta (str): Path to the original FASTA file (CNHEVSAPs_characterized.faa).
        output_fasta (str): Path for the new FASTA file containing only enzyme sequences (CNHEVSACE.faa).
    """
    
    print(f"--- Starting Enzyme Extraction Process (Corrected) ---")
    
    # --- 1. Identify Enzyme Accession IDs from InterProScan results ---
    
    if not os.path.exists(interpro_file):
        print(f"Error: InterProScan file not found at '{interpro_file}'. Please check the path.")
        return

    # CORRECTED: Use usecols=[0, 12]. Column 0 is Accession ID. Column 12 is InterPro Description.
    try:
        df = pd.read_csv(
            interpro_file, 
            sep='\t', 
            header=None, 
            usecols=[0, 12],
            names=['Accession_ID', 'InterPro_Description'],
            # Ensure proper handling of missing values in the description column
            na_filter=False 
        )
        
    except Exception as e:
        print(f"Error reading TSV file: {e}")
        print("Please ensure the file format is correct (tab-separated) and has the expected columns.")
        return

    print(f"Successfully loaded {len(df)} lines from InterProScan results.")

    # A comprehensive list of enzyme classification terms to ensure we catch all types
    enzyme_keywords = [
        'Enzyme', 'EC:', 'transferase', 'lyase', 'hydrolase', 'isomerase', 
        'ligase', 'oxidoreductase', 'synthase', 'kinase', 'phosphatase', 
        'reductase', 'dehydrogenase'
    ]
    
    # Create a regular expression pattern for case-insensitive search
    search_pattern = '|'.join(enzyme_keywords)
    
    # Filter the DataFrame: check if the Description contains any of the enzyme keywords
    is_enzyme = df['InterPro_Description'].astype(str).str.contains(
        search_pattern, 
        case=False
    )
    
    # Get the unique set of accession IDs for the identified enzymes
    enzyme_accessions = set(df[is_enzyme]['Accession_ID'])
    
    print(f"Identified {len(enzyme_accessions)} unique proteins classified as enzymes.")

    # --- 2. Fetch Sequences for Identified Enzyme Accessions ---

    if not os.path.exists(input_fasta):
        # Fallback check for the expected location based on your previous description
        print(f"Warning: Input FASTA file not found at '{input_fasta}'. Checking 'output_faa/{input_fasta}'.")
        input_fasta = os.path.join("output_faa", input_fasta)
        
        if not os.path.exists(input_fasta):
            print(f"Error: Input FASTA file not found at '{input_fasta}'. Please verify the location.")
            return
        
    extracted_sequences = []
    found_count = 0
    
    # Parse the input FASTA file
    for record in SeqIO.parse(input_fasta, "fasta"):
        # We need the accession ID from the FASTA header (before the first space)
        # InterProScan often uses the full header as the ID, so we check both the full ID and the first word.
        # Check if the full ID is present
        if record.id in enzyme_accessions:
            extracted_sequences.append(record)
            found_count += 1
            enzyme_accessions.discard(record.id)
        # Check if only the first word/accession part is present (e.g., if the InterPro ID was trimmed)
        elif record.id.split()[0] in enzyme_accessions:
            extracted_sequences.append(record)
            found_count += 1
            enzyme_accessions.discard(record.id.split()[0])
        

    print(f"Successfully extracted {found_count} enzyme sequences from the FASTA file.")

    # --- 3. Save the Sequences to the new FASTA file ---
    
    with open(output_fasta, "w") as f:
        SeqIO.write(extracted_sequences, f, "fasta")

    print(f"Completed! Enzyme sequences saved to '{output_fasta}'.")

    if enzyme_accessions:
        print(f"\nNote: {len(enzyme_accessions)} enzyme IDs were identified in the TSV but could not be matched to sequences in the FASTA file.")


# --- Example Usage ---
if __name__ == "__main__":
    interpro_file = "InterProScan_result.tsv"
    # Assuming the input FASTA file is CNHEVSAPs_characterized.faa (in the current directory)
    input_fasta_file = "CHNEV_cellular_location.faa" 
    output_fasta_file = "CNHEVCME.faa"
    
    extract_enzyme_targets(interpro_file, input_fasta_file, output_fasta_file)
