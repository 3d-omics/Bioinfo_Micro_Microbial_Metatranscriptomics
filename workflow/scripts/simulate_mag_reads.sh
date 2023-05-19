#!/usr/bin/env bash

set -euo pipefail

sample_ids=(GBRF1.1 GBRF2.1 GBRF3.1 GBRM1.1 GBRM2.1 GBRM3.1)
mag_ids=(acutalibacter_ornithocaccae scybalocola_faecipullorum)


mkdir --parents resources/reads_mags resources/reads

for sample_id in "${sample_ids[@]}" ; do

    for mag_id in "${mag_ids[@]}"; do

        wgsim \
            -N 1000 \
            -1 150 \
            -2 150 \
            resources/mags/"$mag_id".fna.gz \
            >(gzip -9 > resources/reads_mags/"$sample_id"_"$mag_id"_1.fq.gz) \
            >(gzip -9 > resources/reads_mags/"$sample_id"_"$mag_id"_2.fq.gz) \

    done

    cat \
        resources/reads_host/"$sample_id"_sub_1.fq.gz \
        resources/reads_mags/"$sample_id"_*_1.fq.gz \
    > resources/reads_mixed/"$sample_id"_1.fq.gz

    cat \
        resources/reads_host/"$sample_id"_sub_2.fq.gz \
        resources/reads_mags/"$sample_id"_*_2.fq.gz \
    > resources/reads_mixed/"$sample_id"_2.fq.gz

done
