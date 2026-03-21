#!/bin/bash

OUTDIR="bracken_results"
mkdir -p "$OUTDIR"

for REPORT in kraken_results/*.report

do
SAMPLE=$(basename "$REPORT" .report)

bracken \
	-d db/ \
	-i "$REPORT" \
	-o "${OUTDIR}/${SAMPLE}.bracken" \
	-t 0
done
