#!/usr/bin/env bash
set -euo pipefail


tail -n+2 "$@" \
| awk \
    -v OFS="\t" \
    '{print $3, $5 - 1, $6, $1, 1000, ($7 == 1 ? "+" : "-")}' \
| bedtools sort
