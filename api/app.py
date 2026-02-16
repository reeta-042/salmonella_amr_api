from fastapi import FastAPI, UploadFile, File, HTTPException, status
from fastapi.responses import JSONResponse
from fastapi.middleware.cors import CORSMiddleware
import subprocess
import pickle
import pandas as pd
import numpy as np
import shap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import uuid
import shutil
import base64
from io import BytesIO
from datetime import datetime
from pathlib import Path
import traceback

app = FastAPI(
    title="Salmonella AMR Prediction API",
    version="1.0.0",
    description="Predict antibiotic resistance in Salmonella genomes using ML with SHAP explainability"
)

# CORS for frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

WORK_DIR = Path("/app/work")
MODELS_DIR = Path(os.getenv('MODELS_DIR', '/app/models'))
SCRIPTS_DIR = Path(os.getenv('SCRIPTS_DIR', '/app/scripts'))

MODELS = {
    "pefoxacin_full": MODELS_DIR / "pefoxacin_full_model.pkl",
    "trimethoprim_full": MODELS_DIR / "trimethoprim_full_model.pkl",
    "sulfamethoxazole_full": MODELS_DIR / "sulfamethoxazole_full_model.pkl",
    "pefoxacin_partial": MODELS_DIR / "pefoxacin_snps_kmers_model.pkl",
    "trimethoprim_partial": MODELS_DIR / "trimethoprim_snps_kmers_model.pkl",
    "sulfamethoxazole_partial": MODELS_DIR / "sulfamethoxazole_genes_snps_model.pkl"
}

# Global exception handler
@app.exception_handler(Exception)
async def global_exception_handler(request, exc):
    return JSONResponse(
        status_code=500,
        content={
            "status": "error",
            "message": str(exc),
            "type": type(exc).__name__,
            "timestamp": datetime.utcnow().isoformat() + "Z"
        }
    )

def get_shap_explanation(model, X, top_n=5):
    """Extract top contributing features using SHAP"""
    try:
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
            impact = row['shap_value']
            direction = "promotes_resistance" if impact > 0 else "promotes_susceptibility"
            
            feat_name = row['feature']
            if '>' in feat_name and 'NC_' in feat_name:
                feat_type = "SNP"
            elif len(feat_name) == 10 and feat_name.isalpha():
                feat_type = "K-mer"
            else:
                feat_type = "Gene"
            
            evidence.append({
                "feature": feat_name,
                "type": feat_type,
                "impact_score": round(float(impact), 4),
                "effect": direction
            })
        
        return evidence
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"SHAP analysis failed: {str(e)}")

def create_force_plot(model, X, antibiotic):
    """Generate SHAP force plot as base64 image"""
    try:
        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(X)
        
        if isinstance(shap_values, list):
            shap_values = shap_values[1]
            expected_value = explainer.expected_value[1]
        else:
            expected_value = explainer.expected_value
        
        plt.figure(figsize=(14, 3))
        shap.force_plot(
            expected_value,
            shap_values[0],
            X.iloc[0],
            matplotlib=True,
            show=False,
            text_rotation=10
        )
        plt.title(f"{antibiotic} - Feature Contributions", fontsize=14, fontweight='bold')
        
        buf = BytesIO()
        plt.savefig(buf, format='png', dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()
        buf.seek(0)
        
        img_base64 = base64.b64encode(buf.read()).decode('utf-8')
        return f"data:image/png;base64,{img_base64}"
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"SHAP visualization failed: {str(e)}")

@app.get("/health")
def health_check():
    """Health check endpoint"""
    try:
        # Verify models exist
        models_ok = all(p.exists() for p in MODELS.values())
        
        return {
            "status": "healthy" if models_ok else "degraded",
            "version": "1.0.0",
            "models_loaded": models_ok,
            "timestamp": datetime.utcnow().isoformat() + "Z"
        }
    except Exception as e:
        return JSONResponse(
            status_code=503,
            content={
                "status": "unhealthy",
                "error": str(e),
                "timestamp": datetime.utcnow().isoformat() + "Z"
            }
        )

@app.post("/predict")
async def predict_resistance(genome: UploadFile = File(...)):
    """
    Predict antibiotic resistance for Salmonella genome.
    
    Args:
        genome: FASTA/FNA file (assembled genome, 1MB-10MB)
    
    Returns:
        JSON with predictions for pefoxacin, trimethoprim, sulfamethoxazole
    """
    job_id = str(uuid.uuid4())[:8]
    work_path = WORK_DIR / job_id
    
    try:
        # Validate file format
        if not genome.filename.endswith(('.fna', '.fasta', '.fa')):
            raise HTTPException(
                status_code=400,
                detail="Invalid file format. Please upload .fna, .fasta, or .fa file"
            )
        
        # Read and validate file size
        content = await genome.read()
        file_size = len(content)
        
        if file_size > 10 * 1024 * 1024:  # 10MB
            raise HTTPException(
                status_code=400,
                detail="File too large. Maximum size is 10MB"
            )
        
        if file_size < 1000:  # 1KB minimum
            raise HTTPException(
                status_code=400,
                detail="File too small. Please upload a valid genome assembly"
            )
        
        # Create work directory
        work_path.mkdir(parents=True, exist_ok=True)
        start_time = datetime.utcnow()
        
        # Save genome
        genome_path = work_path / "query_genome.fna"
        with open(genome_path, "wb") as f:
            f.write(content)
        
        # Run pipeline
        os.chdir(work_path)
        
        try:
            subprocess.run(["bash", f"{SCRIPTS_DIR}/01_extract_genes.sh"], 
                          check=True, capture_output=True, timeout=300)
            subprocess.run(["python3", f"{SCRIPTS_DIR}/01b_process_genes.py"], 
                          check=True, capture_output=True, timeout=60)
            subprocess.run(["bash", f"{SCRIPTS_DIR}/02_create_blast_db.sh"], 
                          check=True, capture_output=True, timeout=60)
            subprocess.run(["python3", f"{SCRIPTS_DIR}/03_extract_kmers.py"], 
                          check=True, capture_output=True, timeout=600)
            subprocess.run(["bash", f"{SCRIPTS_DIR}/04_extract_snps.sh"], 
                          check=True, capture_output=True, timeout=600)
            subprocess.run(["python3", f"{SCRIPTS_DIR}/04b_process_snps.py"], 
                          check=True, capture_output=True, timeout=60)
            subprocess.run(["python3", f"{SCRIPTS_DIR}/06_align_features.py"], 
                          check=True, capture_output=True, timeout=60)
        except subprocess.TimeoutExpired as e:
            raise HTTPException(
                status_code=504,
                detail=f"Pipeline timeout during {e.cmd[1].split('/')[-1]}"
            )
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr.decode() if e.stderr else "Unknown error"
            raise HTTPException(
                status_code=500,
                detail=f"Pipeline failed: {error_msg[:200]}"
            )
        
        # Load features
        try:
            df_full = pd.read_csv(work_path / "aligned_full.csv")
            df_snps_kmers = pd.read_csv(work_path / "aligned_snps_kmers.csv")
            df_genes_snps = pd.read_csv(work_path / "aligned_genes_snps.csv")
            genes_df = pd.read_csv(work_path / "gene_presence_production.csv")
            kmers_df = pd.read_csv(work_path / "kmer_production.csv")
            snps_df = pd.read_csv(work_path / "snp_production.csv")
        except FileNotFoundError as e:
            raise HTTPException(
                status_code=500,
                detail=f"Feature extraction incomplete: {str(e)}"
            )
        
        genome_id = df_full['Genome_ID'].iloc[0]
        X_full = df_full.drop(columns=['Genome_ID'])
        X_snps_kmers = df_snps_kmers.drop(columns=['Genome_ID'])
        X_genes_snps = df_genes_snps.drop(columns=['Genome_ID'])
        
        # Predict with SHAP
        results = {}
        
        for antibiotic in ["pefoxacin", "trimethoprim", "sulfamethoxazole"]:
            try:
                with open(MODELS[f"{antibiotic}_full"], 'rb') as f:
                    model_full = pickle.load(f)
                with open(MODELS[f"{antibiotic}_partial"], 'rb') as f:
                    model_partial = pickle.load(f)
            except Exception as e:
                raise HTTPException(
                    status_code=500,
                    detail=f"Failed to load {antibiotic} models: {str(e)}"
                )
            
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
                
                if final_prob >= 0.85:
                    confidence = "High"
                    action = "REPORT_FINAL"
                elif final_prob >= 0.65:
                    confidence = "Medium"
                    action = "CONSIDER_CONFIRMATION"
                else:
                    confidence = "Low"
                    action = "CONFIRMATORY_AST_REQUIRED"
                
                consensus = f"Both models agree ({final_prob:.1%} confident)"
            
            evidence = get_shap_explanation(best_model, best_X, top_n=5)
            force_plot = create_force_plot(best_model, best_X, antibiotic.capitalize())
            
            results[antibiotic] = {
                "phenotype": final_pred,
                "probability_score": round(float(final_prob), 4),
                "confidence_category": confidence,
                "action_required": action,
                "consensus": consensus,
                "evidence": evidence,
                "shap_visualization": force_plot,
                "model_breakdown": {
                    "full_model": {"prediction": full_pred, "probability": round(float(proba_full[1]), 4)},
                    "partial_model": {"prediction": partial_pred, "probability": round(float(proba_partial[1]), 4)}
                }
            }
        
        end_time = datetime.utcnow()
        
        output = {
            "job_id": job_id,
            "sample_id": genome_id,
            "status": "completed",
            "timestamps": {
                "submitted_at": start_time.isoformat() + "Z",
                "completed_at": end_time.isoformat() + "Z",
                "processing_time_seconds": int((end_time - start_time).total_seconds())
            },
            "quality_metrics": {
                "genome_size_mb": round(file_size / (1024 * 1024), 2),
                "genes_detected": int(genes_df.shape[1] - 1),
                "kmers_matched": int(kmers_df.shape[1] - 1),
                "snps_detected": int(snps_df.shape[1] - 1)
            },
            "predictions": results,
            "model_metadata": {
                "pipeline_version": "1.0.0",
                "trained_date": "2026-01-15",
                "training_samples": 338,
                "card_version": "2024.01",
                "reference_genome": "Salmonella_Typhimurium_LT2"
            }
        }
        
        # Cleanup
        shutil.rmtree(work_path)
        
        return output
    
    except HTTPException:
        raise
    except Exception as e:
        if work_path.exists():
            shutil.rmtree(work_path)
        raise HTTPException(
            status_code=500,
            detail=f"Unexpected error: {str(e)}\n{traceback.format_exc()[:500]}"
        )

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)