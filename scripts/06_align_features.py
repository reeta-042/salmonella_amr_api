#!/usr/bin/env python3
import os, sys, pandas as pd

WORK_DIR = os.getcwd()
LOG_FILE = os.path.join(WORK_DIR, "logs", "06_align_features.log")
TEMPLATE_DIR = os.getenv('FEATURE_TEMPLATES_DIR', '/app/feature_templates')

def log(msg):
    print(msg)
    with open(LOG_FILE, 'a') as f:
        f.write(msg + '\n')

with open(LOG_FILE, 'w') as f:
    f.write("SCRIPT 6: FEATURE ALIGNMENT\n")

log("Loading extracted features...")
genes_df = pd.read_csv('gene_presence_production.csv')
kmers_df = pd.read_csv('kmer_production.csv')
snps_df = pd.read_csv('snp_production.csv')

log(f"  Genes: {genes_df.shape[1] - 1}, K-mers: {kmers_df.shape[1] - 1}, SNPs: {snps_df.shape[1] - 1}")

genome_id = genes_df['Genome_ID'].iloc[0]
merged = genes_df.copy()

if kmers_df.shape[1] > 1:
    merged = merged.merge(kmers_df, on='Genome_ID', how='left')
if snps_df.shape[1] > 1:
    merged = merged.merge(snps_df, on='Genome_ID', how='left')

merged.fillna(0, inplace=True)
for col in [c for c in merged.columns if c != 'Genome_ID']:
    merged[col] = merged[col].astype(int)

log(f"  Total merged features: {merged.shape[1] - 1}")

def align_to_template(merged_df, template_file, output_file, template_name):
    log(f"\n  Aligning to: {template_name}")
    
    if not os.path.exists(template_file):
        log(f"✗ ERROR: Template not found: {template_file}")
        sys.exit(1)
    
    with open(template_file, 'r') as f:
        template_features = [line.strip() for line in f if line.strip()]
    
    log(f"    Template features: {len(template_features)}")
    
    aligned_row = {'Genome_ID': genome_id}
    found = 0
    
    for feat in template_features:
        if feat in merged_df.columns:
            aligned_row[feat] = int(merged_df[feat].iloc[0])
            found += 1
        else:
            aligned_row[feat] = 0
    
    pd.DataFrame([aligned_row]).to_csv(output_file, index=False)
    
    missing_pct = (1 - found / len(template_features)) * 100
    log(f"    Matched: {found}/{len(template_features)} ({100-missing_pct:.1f}%)")
    log(f"    ✓ Saved: {output_file}")
    return missing_pct

align_to_template(merged, os.path.join(TEMPLATE_DIR, 'features_full_dataset.txt'), 
                  'aligned_full.csv', 'Full Dataset')

align_to_template(merged, os.path.join(TEMPLATE_DIR, 'features_snps_kmers.txt'), 
                  'aligned_snps_kmers.csv', 'SNPs + K-mers')

align_to_template(merged, os.path.join(TEMPLATE_DIR, 'features_genes_snps.txt'), 
                  'aligned_genes_snps.csv', 'Genes + SNPs')

log("\n✓ All alignments complete")