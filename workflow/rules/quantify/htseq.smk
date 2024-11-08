rule quantify__htseq__count:
    input:
        cram=BOWTIE2 / "{mag_catalogue}.{sample_id}.{library_id}.cram",
        crai=BOWTIE2 / "{mag_catalogue}.{sample_id}.{library_id}.cram.crai",
        reference=MAGS / "{mag_catalogue}.fa.gz",
        fai=MAGS / "{mag_catalogue}.fa.gz.fai",
        annotation=MAGS / "{mag_catalogue}.gtf",
    output:
        counts=HTSEQ / "{mag_catalogue}" / "{sample_id}.{library_id}.tsv.gz",
    log:
        HTSEQ / "{mag_catalogue}" / "{sample_id}.{library_id}.log",
    conda:
        "../../environments/htseq.yml"
    params:
        sample_library=lambda w: f"{w.sample_id}.{w.library_id}",
    shell:
        """
        ( echo -E "gene_id\t{params.sample_library}" \
        | gzip \
        > {output.counts} )
        2> {log}

        ( htseq-count \
            --order pos \
            --type CDS \
            --idattr gene_id \
            {input.cram} \
            {input.annotation} \
        | gzip \
        > {output.counts} \
        ) 2>> {log}
        """


rule quantify__htseq__count__aggregate:
    input:
        tsvs=get_tsvs_for_htseq,
    output:
        HTSEQ / "{mag_catalogue}.tsv.gz",
    log:
        HTSEQ / "{mag_catalogue}.log",
    conda:
        "../../environments/htseq.yml"
    params:
        input_folder=lambda w: HTSEQ / f"{w.mag_catalogue}",
    shell:
        """
        Rscript --vanilla workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_folder} \
            --output-file {output} \
        2> {log}
        """


rule quantify__htseq__all:
    input:
        [HTSEQ / f"{mag_catalogue}.tsv.gz" for mag_catalogue in MAG_CATALOGUES],
