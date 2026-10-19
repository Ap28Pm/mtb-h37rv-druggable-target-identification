#!/usr/bin/env python3
import sys
import os
import csv
from collections import defaultdict
from Bio import SeqIO
import pandas as pd

if len(sys.argv) < 4:
    print("Usage: parse_core_blast.py <blast_tab> <query_fasta> <core_fasta>")
    sys.exit(1)

blast_tab = sys.argv[1]
query_fasta = sys.argv[2]
core_fasta = sys.argv[3]

OUTDIR = "blast_results"
os.makedirs(OUTDIR, exist_ok=True)
out_csv = os.path.join(OUTDIR, "query_vs_core_summary.csv")
out_report = os.path.join(OUTDIR, "query_vs_core_report.txt")
out_matches_faa = "output_faa/core_matches.faa"
os.makedirs(os.path.dirname(out_matches_faa), exist_ok=True)

# thresholds (you can change these values before running if desired)
PIDENT_THRESH = 90.0    # percent identity threshold for "highly-matching"
QCOV_THRESH = 0.8       # query coverage threshold (fraction)
EVALUE_MAX = 1e-5

# load sequence lengths for query and core
def load_lengths(fasta):
    d = {}
    for r in SeqIO.parse(fasta, "fasta"):
        # use record.id (part before first space) as key and also full header
        d[r.id] = len(r.seq)
        # also store first token to be safe (if headers have extra text)
        first = r.id.split()[0]
        if first != r.id:
            d[first] = len(r.seq)
    return d

print("Loading sequence lengths...")
query_lengths = load_lengths(query_fasta)
core_lengths = load_lengths(core_fasta)

# parse BLAST tabular file
if not os.path.exists(blast_tab):
    print("ERROR: BLAST tab file not found:", blast_tab)
    sys.exit(1)

hits = []
with open(blast_tab) as fh:
    for line in fh:
        parts = line.strip().split("\t")
        if len(parts) < 12:
            continue
        qseqid, sseqid, pident, aln_len, qlen, slen, qstart, qend, sstart, send, evalue, bitscore = parts[:12]
        pident = float(pident)
        aln_len = int(aln_len)
        evalue = float(evalue)
        qlen = int(qlen) if qlen.isdigit() else query_lengths.get(qseqid, None)
        slen = int(slen) if slen.isdigit() else core_lengths.get(sseqid, None)
        if qlen is None:
            # try splitting id token
            qlen = query_lengths.get(qseqid.split()[0], None)
        if slen is None:
            slen = core_lengths.get(sseqid.split()[0], None)
        qcov = aln_len / qlen if qlen else None
        scov = aln_len / slen if slen else None
        hits.append({
            "qseqid": qseqid, "sseqid": sseqid, "pident": pident,
            "aln_len": aln_len, "q_len": qlen, "s_len": slen,
            "qcov": qcov, "scov": scov,
            "evalue": evalue, "bitscore": float(bitscore)
        })

df = pd.DataFrame(hits)
if df.empty:
    print("No BLAST hits found.")
    df.to_csv(out_csv, index=False)
    sys.exit(0)

# compute best hit per query by bitscore
best_hits = df.sort_values(["qseqid","bitscore"], ascending=[True, False]).groupby("qseqid").first().reset_index()

# annotate whether best hit meets the high-match thresholds
best_hits["high_match"] = ((best_hits["pident"] >= PIDENT_THRESH) & 
                           (best_hits["qcov"] >= QCOV_THRESH) & 
                           (best_hits["evalue"] <= EVALUE_MAX))

# write CSV of all hits and best hits
df.to_csv(out_csv, index=False)

best_out = os.path.join(OUTDIR, "query_vs_core_best_hits.csv")
best_hits.to_csv(best_out, index=False)

# extract matched core sequences (unique sseqid from best_hits that meet at least evalue threshold)
qualified_sseq = set(best_hits.loc[best_hits["evalue"] <= EVALUE_MAX, "sseqid"])
# also optionally require identity/coverage: keep all with evalue but you can filter further
qualified_sseq_high = set(best_hits.loc[best_hits["high_match"], "sseqid"])

# load core fasta into dict for extraction
core_dict = {r.id: r for r in SeqIO.parse(core_fasta, "fasta")}
# try to match by first token as well
core_by_first = {}
for r in SeqIO.parse(core_fasta, "fasta"):
    core_by_first[r.id.split()[0]] = r

extracted = []
for sid in qualified_sseq:
    if sid in core_dict:
        extracted.append(core_dict[sid])
    elif sid.split()[0] in core_by_first:
        extracted.append(core_by_first[sid.split()[0]])
# write all evalue-qualified matches
SeqIO.write(extracted, out_matches_faa, "fasta")

# report
total_queries = len(set(df["qseqid"]))
queries_with_hits = len(set(df["qseqid"]))
queries_with_high = best_hits["high_match"].sum()
unique_core_matches = len({h for h in df["sseqid"]})

with open(out_report, "w") as rep:
    rep.write("Query vs Core BLAST report\n")
    rep.write("=========================\n")
    rep.write(f"BLAST tab: {blast_tab}\n")
    rep.write(f"Query FASTA: {query_fasta}\n")
    rep.write(f"Core FASTA: {core_fasta}\n\n")
    rep.write(f"Total distinct query sequences in BLAST: {total_queries}\n")
    rep.write(f"Distinct query sequences with >=1 hit: {queries_with_hits}\n")
    rep.write(f"Distinct query sequences with BEST-HIT high-match (pident≥{PIDENT_THRESH}, qcov≥{QCOV_THRESH*100:.0f}%, e≤{EVALUE_MAX}): {queries_with_high}\n")
    rep.write(f"Unique core (subject) sequences observed in hits: {unique_core_matches}\n")
    rep.write(f"Core matches extracted to: {out_matches_faa}\n")
    rep.write(f"All hits CSV: {out_csv}\n")
    rep.write(f"Best-hits CSV: {best_out}\n")

# print concise summary to stdout
print("=== Query vs Core BLAST summary ===")
print(f"Total distinct query sequences present in BLAST results: {total_queries}")
print(f"Best-hit HIGH-MATCH count (pident≥{PIDENT_THRESH}, qcov≥{QCOV_THRESH*100:.0f}%, e≤{EVALUE_MAX}): {queries_with_high}")
print(f"Unique core sequences matched (any evalue): {unique_core_matches}")
print(f"Extracted core matches saved to: {out_matches_faa}")
print(f"Full hits CSV: {out_csv}")
print(f"Best-hits CSV: {best_out}")
print(f"Detailed report: {out_report}")

# If user wants to check a particular query header presence:
print("\nTo check whether a specific query is part of the core, search best-hits CSV (column qseqid -> sseqid) or inspect", best_out)
