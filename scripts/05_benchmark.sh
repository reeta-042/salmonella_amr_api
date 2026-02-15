#!/bin/bash
WORK_DIR=$(pwd)
LOG_DIR="$WORK_DIR/logs"
LOG_FILE="$LOG_DIR/05_benchmark_results.log"
RESULTS_CSV="$WORK_DIR/benchmark_results.csv"
REFERENCE=${REFERENCE_GENOME:-"/app/reference/salmonella_LT2.gbff"}
SCRIPTS_DIR=${SCRIPTS_DIR:-"/app/scripts"}

mkdir -p "$LOG_DIR"

echo "============================================================" | tee "$LOG_FILE"
echo "SCRIPT 5: BENCHMARKING" | tee -a "$LOG_FILE"
echo "Started: $(date)" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"

echo "Combination,Gene_Time_s,BLAST_DB_Time_s,Kmer_Time_s,SNP_Time_s,Total_Time_s" > "$RESULTS_CSV"

echo "[TIMING] Gene Extraction..." | tee -a "$LOG_FILE"
GENE_START=$(date +%s)
abricate --db card query_genome.fna > /tmp/bench_card.tsv 2>/dev/null
abricate --db resfinder query_genome.fna > /tmp/bench_res.tsv 2>/dev/null
abricate --summary /tmp/bench_card.tsv /tmp/bench_res.tsv > /tmp/bench_gene_summary.tsv 2>/dev/null
python3 "$SCRIPTS_DIR/01b_process_genes.py" > /dev/null 2>&1
GENE_END=$(date +%s)
GENE_TIME=$((GENE_END - GENE_START))
echo "  Gene extraction: ${GENE_TIME}s" | tee -a "$LOG_FILE"

echo "[TIMING] BLAST Database..." | tee -a "$LOG_FILE"
BLASTDB_START=$(date +%s)
mkdir -p blast_db
makeblastdb -in query_genome.fna -dbtype nucl -out blast_db/query_genome -parse_seqids > /dev/null 2>&1
BLASTDB_END=$(date +%s)
BLASTDB_TIME=$((BLASTDB_END - BLASTDB_START))
echo "  BLAST DB: ${BLASTDB_TIME}s" | tee -a "$LOG_FILE"

echo "[TIMING] K-mer Extraction..." | tee -a "$LOG_FILE"
KMER_START=$(date +%s)
python3 "$SCRIPTS_DIR/03_extract_kmers.py" > /dev/null 2>&1
KMER_END=$(date +%s)
KMER_TIME=$((KMER_END - KMER_START))
echo "  K-mers: ${KMER_TIME}s" | tee -a "$LOG_FILE"

echo "[TIMING] SNP Extraction..." | tee -a "$LOG_FILE"
SNP_START=$(date +%s)
snippy --outdir snippy_bench_out --ref "$REFERENCE" --ctgs query_genome.fna --cpus 4 --force > /dev/null 2>&1
python3 "$SCRIPTS_DIR/04b_process_snps.py" > /dev/null 2>&1
SNP_END=$(date +%s)
SNP_TIME=$((SNP_END - SNP_START))
echo "  SNPs: ${SNP_TIME}s" | tee -a "$LOG_FILE"

echo "" | tee -a "$LOG_FILE"
echo "RESULTS:" | tee -a "$LOG_FILE"
TOTAL_GENES=$GENE_TIME
TOTAL_SNPS=$((GENE_TIME + SNP_TIME))
TOTAL_KMERS=$((GENE_TIME + BLASTDB_TIME + KMER_TIME))
TOTAL_GENES_SNPS=$((GENE_TIME + SNP_TIME))
TOTAL_GENES_KMERS=$((GENE_TIME + BLASTDB_TIME + KMER_TIME))
TOTAL_SNPS_KMERS=$((GENE_TIME + BLASTDB_TIME + KMER_TIME + SNP_TIME))
TOTAL_FULL=$((GENE_TIME + BLASTDB_TIME + KMER_TIME + SNP_TIME))

echo "Genes Only: ${TOTAL_GENES}s" | tee -a "$LOG_FILE"
echo "SNPs Only: ${TOTAL_SNPS}s" | tee -a "$LOG_FILE"
echo "K-mers Only: ${TOTAL_KMERS}s" | tee -a "$LOG_FILE"
echo "Genes + SNPs: ${TOTAL_GENES_SNPS}s" | tee -a "$LOG_FILE"
echo "Genes + K-mers: ${TOTAL_GENES_KMERS}s" | tee -a "$LOG_FILE"
echo "SNPs + K-mers: ${TOTAL_SNPS_KMERS}s" | tee -a "$LOG_FILE"
echo "Full Dataset: ${TOTAL_FULL}s" | tee -a "$LOG_FILE"

echo "Genes Only,$GENE_TIME,0,0,0,$TOTAL_GENES" >> "$RESULTS_CSV"
echo "SNPs Only,$GENE_TIME,0,0,$SNP_TIME,$TOTAL_SNPS" >> "$RESULTS_CSV"
echo "K-mers Only,$GENE_TIME,$BLASTDB_TIME,$KMER_TIME,0,$TOTAL_KMERS" >> "$RESULTS_CSV"
echo "Genes + SNPs,$GENE_TIME,0,0,$SNP_TIME,$TOTAL_GENES_SNPS" >> "$RESULTS_CSV"
echo "Genes + K-mers,$GENE_TIME,$BLASTDB_TIME,$KMER_TIME,0,$TOTAL_GENES_KMERS" >> "$RESULTS_CSV"
echo "SNPs + K-mers,$GENE_TIME,$BLASTDB_TIME,$KMER_TIME,$SNP_TIME,$TOTAL_SNPS_KMERS" >> "$RESULTS_CSV"
echo "Full Dataset,$GENE_TIME,$BLASTDB_TIME,$KMER_TIME,$SNP_TIME,$TOTAL_FULL" >> "$RESULTS_CSV"

echo "" | tee -a "$LOG_FILE"
echo "✓ Results saved to: benchmark_results.csv" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"

rm -rf snippy_bench_out /tmp/bench_*.tsv