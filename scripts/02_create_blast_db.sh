#!/bin/bash
WORK_DIR=$(pwd)
LOG_DIR="$WORK_DIR/logs"
LOG_FILE="$LOG_DIR/02_create_blast_db.log"
DB_DIR="$WORK_DIR/blast_db"
mkdir -p "$LOG_DIR" "$DB_DIR"

echo "============================================================" | tee "$LOG_FILE"
echo "SCRIPT 2: CREATE BLAST DATABASE" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"

if [ ! -f "query_genome.fna" ]; then
    echo "✗ ERROR: query_genome.fna not found" | tee -a "$LOG_FILE"
    exit 1
fi

echo "Creating BLAST database..." | tee -a "$LOG_FILE"
makeblastdb -in query_genome.fna -dbtype nucl \
    -out "$DB_DIR/query_genome" -parse_seqids >> "$LOG_FILE" 2>&1

if [ $? -ne 0 ]; then
    echo "✗ ERROR: makeblastdb failed" | tee -a "$LOG_FILE"
    exit 1
fi

echo "  ✓ BLAST DB created" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"