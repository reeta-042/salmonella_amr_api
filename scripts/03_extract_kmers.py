#!/usr/bin/env python3
import os, sys, subprocess, pandas as pd, time
from Bio import SeqIO
from collections import Counter
import re

WORK_DIR = os.getcwd()
LOG_FILE = os.path.join(WORK_DIR, "logs", "03_extract_kmers.log")
CARD_PROTEIN_FILE = os.getenv('CARD_PROTEIN_FILE', '/app/card_db/card_all_proteins.fasta')
BLAST_DB = os.path.join(WORK_DIR, "blast_db/query_genome")
TEMPLATE_DIR = os.getenv('FEATURE_TEMPLATES_DIR', '/app/feature_templates')

K_SIZE = 10
E_VALUE = 1e-5
MIN_IDENTITY = 80
MIN_LENGTH = 50
MAX_PROTEINS_PER_GENE = 3

def log(msg):
    print(msg)
    with open(LOG_FILE, 'a') as f:
        f.write(msg + '\n')

with open(LOG_FILE, 'w') as f:
    f.write("SCRIPT 3: K-MER EXTRACTION\n")

log("Loading training k-mer feature list...")
training_kmers = set()
for feat_file in ['features_genes_kmers.txt', 'features_snps_kmers.txt', 'features_full_dataset.txt']:
    path = os.path.join(TEMPLATE_DIR, feat_file)
    if os.path.exists(path):
        with open(path, 'r') as f:
            for line in f:
                feat = line.strip()
                if len(feat) == 10 and feat.isalpha() and feat.isupper():
                    training_kmers.add(feat)

log(f"  Training k-mers to look for: {len(training_kmers)}")

log("Loading resistance gene list...")
gene_df = pd.read_csv('gene_presence_production.csv')
resistance_genes = [col for col in gene_df.columns if col != 'Genome_ID']
log(f"  Resistance genes: {len(resistance_genes)}")

log("Filtering CARD proteins...")
def normalize_gene_name(gene):
    normalized = re.sub(r'^(bla|BLA)', '', gene)
    normalized = re.sub(r'_\d+$', '', normalized)
    return normalized.upper()

gene_to_proteins = {gene: [] for gene in resistance_genes}
filtered_proteins = []
matched_genes = set()

for record in SeqIO.parse(CARD_PROTEIN_FILE, "fasta"):
    for gene in resistance_genes:
        gene_norm = normalize_gene_name(gene)
        patterns = [rf'\b{re.escape(gene)}\b', rf'\b{re.escape(gene_norm)}\b']
        if '(' in gene:
            no_parens = gene.replace('(', '').replace(')', '').replace("'", '')
            patterns.append(rf'\b{re.escape(no_parens)}\b')
        for pattern in patterns:
            if re.search(pattern, record.description, re.IGNORECASE):
                if len(gene_to_proteins[gene]) < MAX_PROTEINS_PER_GENE:
                    filtered_proteins.append(record)
                    matched_genes.add(gene)
                    gene_to_proteins[gene].append(record.id)
                break

filtered_fasta = "resistance_proteins_production.faa"
SeqIO.write(filtered_proteins, filtered_fasta, "fasta")
log(f"  Matched {len(matched_genes)}/{len(resistance_genes)} genes")
log(f"  Filtered proteins: {len(filtered_proteins)}")

log("Running tBLASTn...")
blast_output = "blast_production.tsv"
t_start = time.time()

subprocess.run([
    "tblastn", "-query", filtered_fasta, "-db", BLAST_DB,
    "-out", blast_output, "-outfmt", "6 qseqid sseqid pident length qseq sseq",
    "-evalue", str(E_VALUE), "-max_target_seqs", "5", "-num_threads", "4"
], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, check=True)

elapsed = int(time.time() - t_start)
log(f"  BLAST completed in {elapsed} seconds")

if not os.path.exists(blast_output) or os.path.getsize(blast_output) == 0:
    log("⚠ No BLAST hits - all k-mers will be 0")
    pd.DataFrame([{'Genome_ID': 'query_genome'}]).to_csv('kmer_production.csv', index=False)
    sys.exit(0)

blast_df = pd.read_csv(blast_output, sep='\t', names=['query_id', 'subject_id', 'pident', 'length', 'query_seq', 'subject_seq'])
blast_df = blast_df[(blast_df['pident'] >= MIN_IDENTITY) & (blast_df['length'] >= MIN_LENGTH)]
log(f"  Passing BLAST hits: {len(blast_df)}")

log("Extracting k-mers...")
found_kmers = Counter()
for _, hit in blast_df.iterrows():
    protein_seq = hit['subject_seq'].replace('-', '')
    for i in range(len(protein_seq) - K_SIZE + 1):
        kmer = protein_seq[i:i+K_SIZE]
        if '*' not in kmer and 'X' not in kmer:
            found_kmers[kmer] += 1

matched_kmers = {k: v for k, v in found_kmers.items() if k in training_kmers}
log(f"  K-mers matching training: {len(matched_kmers)}")

kmer_row = {'Genome_ID': 'query_genome'}
kmer_row.update(matched_kmers)
pd.DataFrame([kmer_row]).to_csv('kmer_production.csv', index=False)

os.remove(blast_output)
os.remove(filtered_fasta)
log("✓ Saved: kmer_production.csv")