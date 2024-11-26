include: "htseq_functions.smk"


rule quantify__htseq__count:
    input:
        bam=BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.bam",
        bai=BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.bam.bai",
        annotation=MAGS / "{mag_catalogue}.gff",
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
            --type CDS,tRNA,rRNA \
            --idattr ID \
            {input.bam} \
            {input.annotation} \
        | gzip \
        > {output.counts} \
        ) 2>> {log}
        """


rule quantify__htseq__count__aggregate:
    input:
        lambda w: [
            HTSEQ / w.mag_catalogue / f"{sample_id}.{library_id}.tsv.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        HTSEQ / "{mag_catalogue}.tsv.gz",
    log:
        HTSEQ / "{mag_catalogue}.log",
    params:
        subcommand="join",
    wrapper:
        "v5.2.1/utils/csvtk"


rule quantify__htseq__all:
    input:
        [HTSEQ / f"{mag_catalogue}.tsv.gz" for mag_catalogue in MAG_CATALOGUES],
