import pandas as pd

# thresholds
IDENTITY = 30
COVERAGE = 0.7

# load human proteome into dict (for seq lengths)
human = {}
with open("human_proteome.faa") as f:
    seq_id, seq = "", ""
    for line in f:
        line=line.strip()
        if line.startswith(">"):
            if seq_id: human[seq_id] = seq
            seq_id = line[1:].split()[0]
            seq = ""
        else:
            seq += line
    if seq_id: human[seq_id] = seq

# load core proteins into dict (for lengths)
core = {}
with open("output_faa/core_filtered.faa") as f:
    seq_id, seq = "", ""
    for line in f:
        line=line.strip()
        if line.startswith(">"):
            if seq_id: core[seq_id] = seq
            seq_id = line[1:].split()[0]
            seq = ""
        else:
            seq += line
    if seq_id: core[seq_id] = seq

# read blast results
cols = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"]
df = pd.read_csv("blast_results/core_vs_human.tsv", sep="\t", names=cols)

# filter by coverage and identity
filtered = []
for _, row in df.iterrows():
    qid, sid = row["qseqid"], row["sseqid"]
    aln_len = row["length"]
    pident = row["pident"]
    evalue = row["evalue"]
    qcov = aln_len / len(core[qid])
    if pident >= IDENTITY and qcov >= COVERAGE and evalue <= 1e-5:
        filtered.append(qid)

filtered = set(filtered)

# proteins non-homologous to human
non_human_core = {k:v for k,v in core.items() if k not in filtered}

# write new fasta
with open("output_faa/core_non_human.faa","w") as out:
    for k,v in non_human_core.items():
        out.write(f">{k}\n{v}\n")

print(f"Total core proteins: {len(core)}")
print(f"Human homologs removed: {len(filtered)}")
print(f"Non-human core proteins: {len(non_human_core)}")
