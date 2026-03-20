#!/bin/bash

OUTDIR="kraken_results"
mkdir -p "$OUTDIR"

for R1 in trimmed/*_r1_trim.fastq

do
SAMPLE=$(basename "$R1" _r1_trim.fastq)
R2="trimmed/${SAMPLE}_r2_trim.fastq"

k2 classify --db db/ \
	 --threads 8 \
	 --output "${OUTDIR}/${SAMPLE}.kraken" \
	 --report "${OUTDIR}/${SAMPLE}.report" \
	 --confidence 0.2 \
	 --memory-mapping \
	 --use-names \
	 --report-zero-counts \
	 --paired "$R1" "$R2"
done 

