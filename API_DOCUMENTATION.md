# Salmonella AMR Prediction API - Documentation

## Overview
ML-powered API for predicting antibiotic resistance in Salmonella genomes with SHAP explainability.

**Version:** 1.0.0  
**Docker Image:** `ejiaborrita/requence_amr_project:v1.0.0`

---

## Deployment

### Pull Image
```bash
docker pull ejiaborrita/requence_amr_project:v1.0.0
```

### Run Container
```bash
docker run -d -p 8000:8000 --name salmonella-amr ejiaborrita/requence_amr_project:v1.0.0
```

### Health Check
```bash
curl http://localhost:8000/health
```

---

## API Endpoints

### 1. Health Check
**GET** `/health`

Returns API health status and model availability.

**Response (200 OK):**
```json
{
  "status": "healthy",
  "version": "1.0.0",
  "models_loaded": true,
  "timestamp": "2026-02-14T20:00:00Z"
}
```

---

### 2. Predict Resistance
**POST** `/predict`

Predicts antibiotic resistance for 3 antibiotics: Pefoxacin, Trimethoprim, Sulfamethoxazole.

**Request:**
- **Content-Type:** `multipart/form-data`
- **Body:**
  - `genome` (file): FASTA/FNA format, 1KB-10MB

**cURL Example:**
```bash
curl -X POST http://localhost:8000/predict \
  -F "genome=@genome.fna"
```

**Processing Time:** 5-10 minutes

**Response (200 OK):**
```json
{
  "job_id": "abc12345",
  "sample_id": "query_genome",
  "status": "completed",
  "timestamps": {
    "submitted_at": "2026-02-14T20:00:00Z",
    "completed_at": "2026-02-14T20:08:30Z",
    "processing_time_seconds": 510
  },
  "quality_metrics": {
    "genome_size_mb": 4.82,
    "genes_detected": 26,
    "kmers_matched": 1486,
    "snps_detected": 3542
  },
  "predictions": {
    "pefoxacin": {
      "phenotype": "Resistant",
      "probability_score": 0.9234,
      "confidence_category": "High",
      "action_required": "REPORT_FINAL",
      "consensus": "Both models agree (92.3% confident)",
      "evidence": [
        {
          "feature": "gyrA_S83F",
          "type": "SNP",
          "impact_score": 0.4512,
          "effect": "promotes_resistance"
        }
      ],
      "shap_visualization": "data:image/png;base64,...",
      "model_breakdown": {
        "full_model": {"prediction": "Resistant", "probability": 0.9200},
        "partial_model": {"prediction": "Resistant", "probability": 0.9500}
      }
    },
    "trimethoprim": {...},
    "sulfamethoxazole": {...}
  },
  "model_metadata": {
    "pipeline_version": "1.0.0",
    "trained_date": "2026-01-15",
    "training_samples": 338,
    "card_version": "2024.01",
    "reference_genome": "Salmonella_Typhimurium_LT2"
  }
}
```

---

## Response Fields

### Prediction Object

| Field | Type | Description |
|-------|------|-------------|
| `phenotype` | string | "Resistant" or "Susceptible" |
| `probability_score` | float | Confidence probability (0.0-1.0) |
| `confidence_category` | string | "High" (≥85%), "Medium" (65-84%), "Low" (<65%) |
| `action_required` | string | Clinical action: "REPORT_FINAL", "CONSIDER_CONFIRMATION", or "CONFIRMATORY_AST_REQUIRED" |
| `consensus` | string | Agreement status between models |
| `evidence` | array | Top 5 SHAP features influencing prediction |
| `shap_visualization` | string | Base64-encoded PNG force plot |
| `model_breakdown` | object | Individual model predictions for transparency |

### SHAP Visualization
The `shap_visualization` field contains a base64-encoded PNG image showing feature contributions.

**To display in frontend:**
```html
<img src="${prediction.shap_visualization}" alt="SHAP Analysis" />
```

**Color Legend:**
- 🔴 Red: Features promoting resistance
- 🔵 Blue: Features promoting susceptibility

---

## Error Responses

### 400 Bad Request
```json
{
  "status": "error",
  "message": "Invalid file format. Please upload .fna, .fasta, or .fa file",
  "type": "HTTPException",
  "timestamp": "2026-02-14T20:00:00Z"
}
```

**Causes:**
- Invalid file format
- File too large (>10MB)
- File too small (<1KB)

### 500 Internal Server Error
```json
{
  "status": "error",
  "message": "Pipeline failed: SNP extraction error",
  "type": "HTTPException",
  "timestamp": "2026-02-14T20:00:00Z"
}
```

**Causes:**
- Pipeline execution failure
- Model loading error
- SHAP analysis error

### 504 Gateway Timeout
```json
{
  "status": "error",
  "message": "Pipeline timeout during 04_extract_snps.sh",
  "type": "HTTPException",
  "timestamp": "2026-02-14T20:00:00Z"
}
```

---

## System Requirements

**Minimum Resources:**
- CPU: 4 cores
- RAM: 8GB
- Disk: 20GB
- Network: Outbound access for Docker Hub pull

**Rate Limiting:**
- Recommended: 1-2 concurrent requests
- Each prediction is CPU/memory intensive

---

## Model Information

**Antibiotics Covered:**
1. Pefoxacin (Quinolone)
2. Trimethoprim (Folate pathway inhibitor)
3. Sulfamethoxazole (Folate pathway inhibitor)

**Training Data:**
- 338 Salmonella genomes
- Trained: January 2026
- Validation F1: 85-99%

**Feature Types:**
- Resistance genes (ABRicate: CARD + ResFinder)
- 10-mer amino acid k-mers (tBLASTn)
- SNPs vs reference (Snippy)

**Reference Genome:**
- Salmonella enterica serovar Typhimurium LT2
- GenBank: NC_003197

---

## Integration Example (Python)

```python
import requests

# Upload genome
with open('genome.fna', 'rb') as f:
    response = requests.post(
        'http://localhost:8000/predict',
        files={'genome': f}
    )

result = response.json()

# Check prediction
for antibiotic, prediction in result['predictions'].items():
    print(f"{antibiotic}: {prediction['phenotype']}")
    print(f"Confidence: {prediction['confidence_category']}")
    print(f"Action: {prediction['action_required']}")
```

---

## Support & Maintenance

**Version Updates:**
- Semantic versioning: `v1.0.0`, `v1.1.0`, `v2.0.0`
- Check Docker Hub for latest: https://hub.docker.com/r/ejiaborrita/requence_amr_project

**Retraining Schedule:**
- CARD database updates: Quarterly
- Model retraining: After 50-100 flagged cases
- Annual maintenance releases

---

## License
Proprietary - Requence Project

## Contact
GitHub: reeta-042  
Docker Hub: ejiaborrita