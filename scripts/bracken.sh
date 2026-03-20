#!/bin/bash

OUTDIR="bracken_unfilt_results"
mkdir -p "$OUTDIR"

for REPORT in kraken_results/*.report

do
SAMPLE=$(basename "$REPORT" .report)

bracken \
	-d db/ \
	-i "$REPORT" \
	-o "${OUTDIR}/${SAMPLE}_unfilt.bracken" \
	-t 0
done
