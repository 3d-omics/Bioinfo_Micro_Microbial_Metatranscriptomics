rule quantify__subread__feature_counts:
    input:
        cram=BOWTIE2 / "{mag_catalogue}.{sample_id}.{library_id}.cram",
        crai=BOWTIE2 / "{mag_catalogue}.{sample_id}.{library_id}.cram.crai",
        reference=MAGS / "{mag_catalogue}.fa.gz",
        fai=MAGS / "{mag_catalogue}.fa.gz.fai",
        annotation=MAGS / "{mag_catalogue}.gtf",
    output:
        counts=SUBREAD / "{mag_catalogue}" / "{sample_id}.{library_id}.tsv.gz",
        summary=temp(
            SUBREAD / "{mag_catalogue}" / "{sample_id}.{library_id}.tsv.summary"
        ),
    log:
        SUBREAD / "{mag_catalogue}" / "{sample_id}.{library_id}.log",
    conda:
        "../../environments/subread.yml"
    params:
        sample_library=lambda w: f"{w.sample_id}.{w.library_id}",
        tmp_out=lambda w: SUBREAD
        / w.mag_catalogue
        / f"{w.sample_id}.{w.library_id}.tsv",
    shell:
        """
        ( samtools view \
            --reference {input.reference} \
            {input.cram} \
        | featureCounts \
            -F GTF \
            -t CDS \
            -g gene_id \
            -p \
            -a {input.annotation} \
            -o {params.tmp_out} \
        ) 2> {log} 1>&2

        ( grep -v ^# {params.tmp_out} \
        | cut -f 1,7 \
        | gzip \
        > {output.counts} \
        ) 2>> {log}

        rm -rfv {params.tmp_out} 2>> {log} 1>&2
        """


rule quantify__subread__aggregate:
    input:
        tsvs=get_tsvs_for_subread,
    output:
        SUBREAD / "{mag_catalogue}.tsv.gz",
    log:
        SUBREAD / "{mag_catalogue}.log",
    conda:
        "../../environments/subread.yml"
    params:
        input_folder=lambda w: SUBREAD / f"{w.mag_catalogue}",
    shell:
        """
        Rscript --vanilla workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_folder} \
            --output-file {output} \
        2> {log}
        """


rule quantify__subread__all:
    input:
        [SUBREAD / f"{mag_catalogue}.tsv.gz" for mag_catalogue in MAG_CATALOGUES],
