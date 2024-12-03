#!/usr/bin/env bash
set -euo pipefail

# Assume kraken2 is installed

# The file should be a fasta (not gzipped) in resources/kraken_mock.fa

mkdir --parents resources/databases/kraken2/

# Build db
gzip -dc resources/reference/mags_mock.fa.gz \
| awk \
    '{ if($1 ~ "mag01") { print $0 "|kraken:taxid|2838416" } else { print $0 } }' \
| awk \
    '{ if($1 ~ "mag02") { print $0 "|kraken:taxid|2840942" } else { print $0 } }' \
> resources/databases/kraken2/kraken2_mock.fa


mkdir --parents resources/databases/kraken2/kraken2_mock

# Download and format db (this is +40GB once decompressed)
kraken2-build \
    --download-taxonomy \
    --db resources/databases/kraken2/kraken2_mock


# Add fasta, not gzipped
kraken2-build \
    `#--no-masking` \
    --add-to-library resources/databases/kraken2/kraken2_mock.fa \
    --db resources/databases/kraken2/kraken2_mock \
    --threads 4

# Build the db itself
kraken2-build \
    --build \
    --db resources/databases/kraken2/kraken2_mock

# Add bracken files
bracken-build \
    -d resources/databases/kraken2/kraken2_mock \
    -t 8 \
    -k 35 \
    -l 100 \
    -y kraken2


# back up just in case
tar cvf - resources/databases/kraken2/kraken2_mock | pigz -1 > kraken2_mock.tar.gz

# Remove crap
kraken2-build \
    --clean \
    --db resources/databases/kraken2/kraken2_mock


# Test with itself
kraken2 \
    --db resources/databases/kraken2/kraken2_mock \
    --threads 4 \
    --report kraken_mock.report \
    --output kraken_mock.out \
    resources/databases/kraken2/kraken2_mock.fa


# Run on mock reads
kraken2 \
    --db resources/databases/kraken2/kraken2_mock \
    --threads 8 \
    --report GBRF1.1.report \
    --output GBRF1.1.out \
    --gzip-compressed \
    --paired \
    resources/reads/GBRF1.1_1.fq.gz \
    resources/reads/GBRF1.1_2.fq.gz


# clean up
rm \
    --force \
    --verbose \
    GBRF1.1.out \
    GBRF1.1.report \
    kraken2_mock.tar.gz \
    resources/databases/kraken2/kraken2_mock.fa \
    kraken_mock.out \
    kraken_mock.report
