#!/usr/bin/env python3
import os, sys, pickle, json, pandas as pd
import shap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime

WORK_DIR = os.getcwd()
MODELS_DIR = os.getenv('MODELS_DIR', '/app/models')
LOG_DIR = os.path.join(WORK_DIR, "logs")
if not os.path.exists(LOG_DIR):
    os.makedirs(LOG_DIR)
LOG_FILE = os.path.join(LOG_DIR, "07_predict.log")

def log(msg):
    timestamp = datetime.now().strftime("%H:%M:%S")
    line = f"[{timestamp}] {msg}"
    print(line)
    with open(LOG_FILE, 'a') as f:
        f.write(line + '\n')

with open(LOG_FILE, 'w') as f:
    f.write("="*60 + '\n')
    f.write("SCRIPT 7: RESISTANCE PREDICTION WITH SHAP\n")
    f.write(f"Started: {datetime.now()}\n")
    f.write("="*60 + '\n')

def get_shap_explanation(model, X, top_n=5):
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X)
    if isinstance(shap_values, list):
        shap_values = shap_values[1]
    
    feature_impacts = pd.DataFrame({
        'feature': X.columns,
        'shap_value': shap_values[0]
    }).sort_values('shap_value', key=abs, ascending=False)
    
    evidence = []
    for _, row in feature_impacts.head(top_n).iterrows():
        feat = row['feature']
        impact = row['shap_value']
        direction = "promotes_resistance" if impact > 0 else "promotes_susceptibility"
        feat_type = "SNP" if '>' in feat and 'NC_' in feat else ("K-mer" if len(feat)==10 and feat.isalpha() else "Gene")
        evidence.append({
            "feature": feat,
            "type": feat_type,
            "impact_score": round(float(impact), 4),
            "effect": direction
        })
    return evidence

log("Validating files...")
required_files = ['aligned_full.csv', 'aligned_snps_kmers.csv', 'aligned_genes_snps.csv']
for f in required_files:
    if not os.path.exists(f):
        log(f"✗ ERROR: {f} not found")
        sys.exit(1)

log("Loading features...")
df_full = pd.read_csv('aligned_full.csv')
df_snps_kmers = pd.read_csv('aligned_snps_kmers.csv')
df_genes_snps = pd.read_csv('aligned_genes_snps.csv')
genes_df = pd.read_csv('gene_presence_production.csv')
kmers_df = pd.read_csv('kmer_production.csv')
snps_df = pd.read_csv('snp_production.csv')

genome_id = df_full['Genome_ID'].iloc[0]
X_full = df_full.drop(columns=['Genome_ID'])
X_snps_kmers = df_snps_kmers.drop(columns=['Genome_ID'])
X_genes_snps = df_genes_snps.drop(columns=['Genome_ID'])

log("Running predictions with SHAP...")
results = {}

for antibiotic in ["pefoxacin", "trimethoprim", "sulfamethoxazole"]:
    with open(os.path.join(MODELS_DIR, f'{antibiotic}_full_model.pkl'), 'rb') as f:
        model_full = pickle.load(f)
    with open(os.path.join(MODELS_DIR, f'{antibiotic}_snps_kmers_model.pkl' if antibiotic != "sulfamethoxazole" else f'{antibiotic}_genes_snps_model.pkl'), 'rb') as f:
        model_partial = pickle.load(f)
    
    X_partial = X_snps_kmers if antibiotic in ["pefoxacin", "trimethoprim"] else X_genes_snps
    
    proba_full = model_full.predict_proba(X_full)[0]
    proba_partial = model_partial.predict_proba(X_partial)[0]
    
    full_pred = "Resistant" if proba_full[1] >= 0.5 else "Susceptible"
    partial_pred = "Resistant" if proba_partial[1] >= 0.5 else "Susceptible"
    
    if full_pred != partial_pred:
        final_pred = full_pred if proba_full[1] > proba_partial[1] else partial_pred
        final_prob = max(proba_full[1], proba_partial[1]) if final_pred == "Resistant" else max(proba_full[0], proba_partial[0])
        confidence = "Low"
        action = "CONFIRMATORY_AST_REQUIRED"
        consensus = f"Models disagree: Full={full_pred} ({proba_full[1]:.1%}), Partial={partial_pred} ({proba_partial[1]:.1%})"
        best_model = model_full if proba_full[1] > proba_partial[1] else model_partial
        best_X = X_full if proba_full[1] > proba_partial[1] else X_partial
    else:
        final_pred = full_pred
        if proba_full[1] > proba_partial[1]:
            final_prob = proba_full[1] if final_pred == "Resistant" else proba_full[0]
            best_model = model_full
            best_X = X_full
        else:
            final_prob = proba_partial[1] if final_pred == "Resistant" else proba_partial[0]
            best_model = model_partial
            best_X = X_partial
        
        confidence = "High" if final_prob >= 0.85 else ("Medium" if final_prob >= 0.65 else "Low")
        action = "REPORT_FINAL" if final_prob >= 0.85 else ("CONSIDER_CONFIRMATION" if final_prob >= 0.65 else "CONFIRMATORY_AST_REQUIRED")
        consensus = f"Both models agree ({final_prob:.1%} confident)"
    
    evidence = get_shap_explanation(best_model, best_X, top_n=5)
    
    log(f"\n{antibiotic.upper()}:")
    log(f"  Phenotype: {final_pred}")
    log(f"  Probability: {final_prob:.1%}")
    log(f"  Confidence: {confidence}")
    log(f"  Action: {action}")
    log(f"  Evidence (top 5 features):")
    for e in evidence:
        log(f"    - {e['feature']} ({e['type']}): {e['impact_score']:.4f} ({e['effect']})")
    
    results[antibiotic] = {
        "phenotype": final_pred,
        "probability_score": round(float(final_prob), 4),
        "confidence_category": confidence,
        "action_required": action,
        "consensus": consensus,
        "evidence": evidence,
        "model_breakdown": {
            "full_model": {"prediction": full_pred, "probability": round(float(proba_full[1]), 4)},
            "partial_model": {"prediction": partial_pred, "probability": round(float(proba_partial[1]), 4)}
        }
    }

output = {
    "genome_id": genome_id,
    "prediction_timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    "quality_metrics": {
        "genes_detected": int(genes_df.shape[1] - 1),
        "kmers_matched": int(kmers_df.shape[1] - 1),
        "snps_detected": int(snps_df.shape[1] - 1)
    },
    "predictions": results,
    "model_metadata": {
        "pipeline_version": "v1.0.0",
        "trained_date": "2026-01-15",
        "training_samples": 338,
        "card_version": "2024.01",
        "reference_genome": "Salmonella_Typhimurium_LT2"
    }
}

with open('prediction_results.json', 'w') as f:
    json.dump(output, f, indent=4)

log("\n✓ Results saved to: prediction_results.json")
log("="*60)
log(f"COMPLETE: {datetime.now()}")
log("="*60)