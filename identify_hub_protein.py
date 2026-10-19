#!/usr/bin/env python3
import pandas as pd

# Input / Output
input_tsv = "string_interactions_short.tsv"
output_csv = "string_hub_proteins.csv"
min_score = 0.7  # Minimum combined_score to consider an interaction

# Read STRING tsv and fix header
df = pd.read_csv(input_tsv, sep="\t", comment=None)
df.columns = [col.lstrip('#') for col in df.columns]  # remove leading '#' from header if present

# Filter by combined_score threshold
df_filtered = df[df['combined_score'] >= min_score]

# Calculate degree (number of interactions per protein)
degree_dict = {}
for node in df_filtered['node1']:
    degree_dict[node] = degree_dict.get(node, 0) + 1
for node in df_filtered['node2']:
    degree_dict[node] = degree_dict.get(node, 0) + 1

# Convert to DataFrame and sort by degree
degree_df = pd.DataFrame(list(degree_dict.items()), columns=['Protein', 'Degree'])
degree_df = degree_df.sort_values(by='Degree', ascending=False)

# Save full table
degree_df.to_csv(output_csv, index=False)

# Print top 10 hub proteins
print("Top 10 hub proteins (high-confidence interactions):")
print(degree_df.head(10))
print(f"\nFull hub protein table saved as: {output_csv}")
