#!/usr/bin/env python3
import os, sys, pandas as pd

WORK_DIR = os.getcwd()
LOG_FILE = os.path.join(WORK_DIR, "logs", "04b_process_snps.log")
SNIPPY_TAB = os.path.join(WORK_DIR, "snippy_production_out", "snps.tab")
TEMPLATE_DIR = os.getenv('FEATURE_TEMPLATES_DIR', '/app/feature_templates')

def log(msg):
    print(msg)
    with open(LOG_FILE, 'a') as f:
        f.write(msg + '\n')

with open(LOG_FILE, 'w') as f:
    f.write("SCRIPT 4b: PROCESS SNPs\n")

if not os.path.exists(SNIPPY_TAB):
    log(f"✗ ERROR: {SNIPPY_TAB} not found")
    sys.exit(1)

log("Loading training SNP feature list...")
training_snps = set()
for feat_file in ['features_snps_only.txt', 'features_genes_snps.txt', 'features_snps_kmers.txt', 'features_full_dataset.txt']:
    path = os.path.join(TEMPLATE_DIR, feat_file)
    if os.path.exists(path):
        with open(path, 'r') as f:
            for line in f:
                feat = line.strip()
                if '>' in feat and feat.startswith('NC_'):
                    training_snps.add(feat)

log(f"  Training SNPs: {len(training_snps)}")

log("Parsing Snippy output...")
snp_tab = pd.read_csv(SNIPPY_TAB, sep='\t', dtype=str)
log(f"  Raw variants: {len(snp_tab)}")

snp_row = {'Genome_ID': 'query_genome'}
matched_snps = 0

for _, row in snp_tab.iterrows():
    try:
        chrom = str(row['CHROM']).strip().split('.')[0]
        pos = str(row['POS']).strip()
        ref = str(row['REF']).strip()
        alt = str(row['ALT']).strip()
        
        if any(v in ['nan', 'NaN', '', '.'] for v in [chrom, pos, ref, alt]):
            continue
        
        feature_name = f"{chrom}_{pos}_{ref}>{alt}"
        if feature_name in training_snps:
            snp_row[feature_name] = 1
            matched_snps += 1
    except:
        continue

log(f"  SNPs matching training: {matched_snps}")
pd.DataFrame([snp_row]).to_csv('snp_production.csv', index=False)
log("✓ Saved: snp_production.csv")