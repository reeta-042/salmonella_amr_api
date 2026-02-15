#!/usr/bin/env python3
import pandas as pd
import os
import sys
from datetime import datetime

WORK_DIR = os.getcwd()
LOG_FILE = os.path.join(WORK_DIR, "logs", "01b_process_genes.log")

def log(msg):
    print(msg)
    with open(LOG_FILE, 'a') as f:
        f.write(msg + '\n')

with open(LOG_FILE, 'w') as f:
    f.write("SCRIPT 1b: PROCESS GENES\n")

log("Loading ABRicate outputs...")

if os.path.exists('query_genome.fna'):
    GENOME_ID = 'query_genome'
else:
    log("✗ ERROR: query_genome.fna not found")
    sys.exit(1)

try:
    summary = pd.read_csv('gene_summary_production.tsv', sep='\t')
    log(f"  Summary matrix shape: {summary.shape}")
except Exception as e:
    log(f"✗ ERROR: {e}")
    sys.exit(1)

log("Processing gene matrix...")

summary.rename(columns={'#FILE': 'Genome_ID'}, inplace=True)
summary['Genome_ID'] = GENOME_ID

def make_binary(val):
    val_str = str(val).strip()
    if val_str in ['0', '', '0.0', '.', 'nan']:
        return 0
    return 1

drop_cols = [c for c in summary.columns if 'NUM_FOUND' in c.upper()]
if drop_cols:
    summary.drop(columns=drop_cols, inplace=True)

skip_cols = ['Genome_ID']
cols_to_binarize = [c for c in summary.columns if c not in skip_cols]
for col in cols_to_binarize:
    summary[col] = summary[col].apply(make_binary)

if len(summary) > 1:
    feature_cols = [c for c in summary.columns if c != 'Genome_ID']
    collapsed = summary[feature_cols].max(axis=0)
    final_row = pd.DataFrame([collapsed])
    final_row.insert(0, 'Genome_ID', GENOME_ID)
    summary = final_row

summary.to_csv('gene_presence_production.csv', index=False)

log(f"  Genome ID: {summary['Genome_ID'].iloc[0]}")
log(f"  Total gene features: {summary.shape[1] - 1}")
log("✓ Saved: gene_presence_production.csv")