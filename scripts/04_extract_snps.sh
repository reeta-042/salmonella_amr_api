#!/bin/bash
WORK_DIR=$(pwd)
LOG_DIR="$WORK_DIR/logs"
LOG_FILE="$LOG_DIR/04_extract_snps.log"
SNIPPY_OUTDIR="$WORK_DIR/snippy_production_out"
REFERENCE="${REFERENCE_GENOME:-/app/reference/salmonella_LT2.gbff}"
mkdir -p "$LOG_DIR"

echo "============================================================" | tee "$LOG_FILE"
echo "SCRIPT 4: SNP EXTRACTION (Snippy)" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"

if [ ! -f "query_genome.fna" ]; then
    echo "✗ ERROR: query_genome.fna not found" | tee -a "$LOG_FILE"
    exit 1
fi

if [ ! -f "$REFERENCE" ]; then
    echo "✗ ERROR: Reference genome not found: $REFERENCE" | tee -a "$LOG_FILE"
    exit 1
fi

echo "✓ Using reference: $REFERENCE" | tee -a "$LOG_FILE"

if [ -d "$SNIPPY_OUTDIR" ]; then
    rm -rf "$SNIPPY_OUTDIR"
fi

echo "[1/2] Running Snippy..." | tee -a "$LOG_FILE"
snippy --outdir "$SNIPPY_OUTDIR" --ref "$REFERENCE" \
    --ctgs query_genome.fna --cpus 4 --force >> "$LOG_FILE" 2>&1

if [ $? -ne 0 ]; then
    echo "✗ ERROR: Snippy failed" | tee -a "$LOG_FILE"
    exit 1
fi

echo "  ✓ Snippy done" | tee -a "$LOG_FILE"

if [ ! -f "$SNIPPY_OUTDIR/snps.tab" ]; then
    echo "✗ ERROR: snps.tab not found" | tee -a "$LOG_FILE"
    exit 1
fi

SNP_COUNT=$(grep -v "^#" "$SNIPPY_OUTDIR/snps.tab" | wc -l)
echo "  ✓ Raw SNPs found: $SNP_COUNT" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"