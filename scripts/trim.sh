#!/bin/bash

OUTDIR="trimmed"
mkdir -p "$OUTDIR"

for R1 in genomes/*_r1.fastq

do
SAMPLE=$(basename "$R1" _r1.fastq)
R2="genomes/${SAMPLE}_r2.fastq"

fastp \
	-i "$R1"\
	-I "$R2" \
	-o "${OUTDIR}/${SAMPLE}_r1_trim.fastq" \
	-O "${OUTDIR}/${SAMPLE}_r2_trim.fastq" \
	-q 20 \
	-h "${OUTDIR}/${SAMPLE}.html" \
	-w 8 && \
	rm "$R1" "$R2"
done

