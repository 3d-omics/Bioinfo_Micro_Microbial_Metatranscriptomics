#!/usr/bin/env bash

set -euo pipefail

export GTDBTK_DATA_PATH="/home/share/db/gtdbtk/release214"


gtdbtk classify_wf \
    --genome_dir separated \
    --extension  fa.gz \
    --out_dir    gtdbtk \
    --cpus       24 \
    --skip_ani_screen






DRAMPATH=/home/share/db/dram/20230811



DRAM-setup.py set_database_locations \
    --amg_database_loc          $DRAMPATH/amg_database.*.tsv \
    --dbcan_fam_activities_loc  $DRAMPATH/CAZyDB.*.fam-activities.txt \
    --dbcan_loc                 $DRAMPATH/dbCAN-HMMdb-V*.txt \
    --dbcan_subfam_ec_loc       $DRAMPATH/CAZyDB.*.fam.subfam.ec.txt \
    --description_db_loc        $DRAMPATH/description_db.sqlite \
    --etc_module_database_loc   $DRAMPATH/etc_mdoule_database.*.tsv \
    --function_heatmap_form_loc $DRAMPATH/function_heatmap_form.*.tsv \
    --genome_summary_form_loc   $DRAMPATH/genome_summary_form.*.tsv \
    --kofam_hmm_loc             $DRAMPATH/kofam_profiles.hmm \
    --kofam_ko_list_loc         $DRAMPATH/kofam_ko_list.tsv \
    --module_step_form_loc      $DRAMPATH/module_step_form.*.tsv \
    --peptidase_loc             $DRAMPATH/peptidases.*.mmsdb \
    --pfam_hmm_loc              $DRAMPATH/Pfam-A.hmm.dat.gz \
    --pfam_loc                  $DRAMPATH/pfam.mmspro \
    --viral_loc                 $DRAMPATH/refseq_viral.*.mmsdb \
    --vog_annotations_loc       $DRAMPATH/vog_annotations_latest.tsv.gz \
    --vogdb_loc                 $DRAMPATH/vog_latest_hmms.txt


mkdir --parents dram/annotate

parallel \
    --jobs 24 \
    --retries 5 \
    DRAM.py annotate \
        --input_fasta {} \
        --output_dir "dram/annotate/{/.}" \
        --threads 1 \
        --gtdb_taxonomy gtdbtk/gtdbtk.bac120.summary.tsv \
::: separated/*.fa.gz

for file in annotations trnas rrnas ; do
    csvstack \
        --tabs \
        dram/annotate/*/*${file}.tsv \
    | csvformat \
        --out-tabs \
    > dram/${file}.tsv
done

DRAM.py distill \
    --input_file dram/annotations.tsv \
    --output_dir dram/distill \
2> dram/distill.log
