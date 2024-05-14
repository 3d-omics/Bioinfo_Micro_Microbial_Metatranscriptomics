#!/usr/bin/env bash
set -euo pipefail


tail -n+2 "$@" \
| awk \
    -v OFS="\t" \
    '{print $3, $5 - 1, $6, $1, 1000, ($7 == 1 ? "+" : "-")}' \
| bedtools sort




bedtools intersect \
    -a test.bam \
    -b annotations.tsv \
    -wb


samtools view -h test.bam | bam2bed - | bedmap --echo --count annotations.bed6


bedtools intersect \
    -a annotations.bed6 \
    -b test.bam \
    -c -bed \
| cut -f 4,7 \
| sort \
| uniq --count \
| awk -v OFS="\t" '{print $2, $1}'
