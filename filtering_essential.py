import pandas as pd

# thresholds (stricter for essential genes)
IDENTITY = 75
COVERAGE = 0.7
EVALUE = 1e-5

# load core_non_human_filtered proteins
core = {}
with open("output_faa/core_non_human.faa") as f:
    seq_id, seq = "", ""
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if seq_id: core[seq_id] = seq
            seq_id = line[1:].split()[0]
            seq = ""
        else:
            seq += line
    if seq_id: core[seq_id] = seq

# read blast results
cols = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"]
df = pd.read_csv("blast_results/DEG_vs_CNH.tsv", sep="\t", names=cols)

# filter hits
essential_hits = set()
for _, row in df.iterrows():
    qid, sid = row["qseqid"], row["sseqid"]
    aln_len = row["length"]
    pident = row["pident"]
    evalue = row["evalue"]
    qcov = aln_len / len(core[qid])
    if pident >= IDENTITY and qcov >= COVERAGE and evalue <= EVALUE:
        essential_hits.add(qid)

# extract essential proteins from core
essential_core = {k:v for k,v in core.items() if k in essential_hits}

# write essential core proteins
with open("CNHE.faa","w") as out:
    for k,v in essential_core.items():
        out.write(f">{k}\n{v}\n")

print(f"Total core non-human proteins: {len(core)}")
print(f"Essential core proteins found: {len(essential_core)}")
