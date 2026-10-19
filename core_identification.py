import glob
import pandas as pd
from collections import defaultdict, Counter

# Thresholds
IDENTITY_THRESH = 90
COVERAGE_THRESH = 0.8

# Load proteomes
proteomes = glob.glob("proteomes/*.faa")
proteome_dict = {}
for f in proteomes:
    name = f.split("/")[-1].replace(".faa","")
    seqs = {}
    with open(f) as handle:
        seq_id, seq = "", ""
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if seq_id: seqs[seq_id] = seq
                seq_id = line[1:].split()[0]
                seq = ""
            else:
                seq += line
        if seq_id: seqs[seq_id] = seq
    proteome_dict[name] = seqs

# Initialize clusters
clusters = []
protein_to_cluster = {}

# Parse BLAST results
blast_files = glob.glob("blast_results/*.tab")
for bf in blast_files:
    with open(bf) as f:
        for line in f:
            cols = line.strip().split("\t")
            qid, sid, pid, aln_len, *_ = cols
            pid = float(pid)
            aln_len = int(aln_len)
            q_len = len(proteome_dict[bf.split('/')[-1].split('_vs_')[0]][qid])
            s_len = len(proteome_dict[bf.split('/')[-1].split('_vs_')[1].replace('.tab','')][sid])
            cov_q = aln_len/q_len
            cov_s = aln_len/s_len
            if pid >= IDENTITY_THRESH and cov_q >= COVERAGE_THRESH and cov_s >= COVERAGE_THRESH:
                # Merge clusters if already exists
                c1 = protein_to_cluster.get(qid)
                c2 = protein_to_cluster.get(sid)
                if c1 and c2:
                    if c1 != c2:
                        clusters[c1].update(clusters[c2])
                        for p in clusters[c2]:
                            protein_to_cluster[p] = c1
                        clusters[c2] = set()
                elif c1:
                    clusters[c1].add(sid)
                    protein_to_cluster[sid] = c1
                elif c2:
                    clusters[c2].add(qid)
                    protein_to_cluster[qid] = c2
                else:
                    idx = len(clusters)
                    clusters.append(set([qid,sid]))
                    protein_to_cluster[qid] = idx
                    protein_to_cluster[sid] = idx

# Add singletons
for genome, seqs in proteome_dict.items():
    for prot in seqs:
        if prot not in protein_to_cluster:
            idx = len(clusters)
            clusters.append(set([prot]))
            protein_to_cluster[prot] = idx

# Build presence/absence matrix
matrix = pd.DataFrame(0, index=range(len(clusters)), columns=proteome_dict.keys())
for i, cluster in enumerate(clusters):
    genomes_in_cluster = set()
    for prot in cluster:
        for g, seqs in proteome_dict.items():
            if prot in seqs:
                matrix.loc[i, g] = 1

matrix.to_csv("core_accessory/presence_absence_matrix.csv", index_label="ClusterID")

# Identify core, accessory, singleton
core = matrix[matrix.sum(axis=1)==len(proteome_dict)]
singleton = matrix[matrix.sum(axis=1)==1]
accessory = matrix[(matrix.sum(axis=1)>1) & (matrix.sum(axis=1)<len(proteome_dict))]

core.to_csv("core_accessory/core_genes.csv")
singleton.to_csv("core_accessory/singleton_genes.csv")
accessory.to_csv("core_accessory/accessory_genes.csv")

# Generate .faa files
def write_faa(cluster_ids, filename):
    with open(filename, "w") as f:
        for cid in cluster_ids:
            for prot in clusters[cid]:
                for g, seqs in proteome_dict.items():
                    if prot in seqs:
                        f.write(f">{g}_{prot}\n{seqs[prot]}\n")

write_faa(core.index, "output_faa/core.faa")
write_faa(accessory.index, "output_faa/accessory.faa")
write_faa(singleton.index, "output_faa/singleton.faa")

# Statistics
report = f"""
Total clusters: {len(clusters)}
Core proteins: {len(core)}
Accessory proteins: {len(accessory)}
Singleton proteins: {len(singleton)}
"""
with open("core_accessory/statistics_report.txt","w") as f:
    f.write(report)
print(report)
