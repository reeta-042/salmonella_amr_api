#!/bin/bash
WORK_DIR=$(pwd)
LOG_DIR="$WORK_DIR/logs"
LOG_FILE="$LOG_DIR/01_extract_genes.log"
mkdir -p "$LOG_DIR"

echo "============================================================" | tee "$LOG_FILE"
echo "SCRIPT 1: GENE EXTRACTION" | tee -a "$LOG_FILE"
echo "Started: $(date)" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"

if [ ! -f "query_genome.fna" ]; then
    echo "✗ ERROR: query_genome.fna not found" | tee -a "$LOG_FILE"
    exit 1
fi

echo "✓ Input genome found" | tee -a "$LOG_FILE"
echo "[1/3] Running ABRicate (CARD)..." | tee -a "$LOG_FILE"
abricate --db card query_genome.fna > card_production.tsv 2>> "$LOG_FILE"
echo "  ✓ CARD done" | tee -a "$LOG_FILE"

echo "[2/3] Running ABRicate (ResFinder)..." | tee -a "$LOG_FILE"
abricate --db resfinder query_genome.fna > resfinder_production.tsv 2>> "$LOG_FILE"
echo "  ✓ ResFinder done" | tee -a "$LOG_FILE"

echo "[3/3] Creating summary..." | tee -a "$LOG_FILE"
abricate --summary card_production.tsv resfinder_production.tsv > gene_summary_production.tsv 2>> "$LOG_FILE"
echo "  ✓ Summary created" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"