rule quantify__subread__feature_counts:
    input:
        bam=BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.bam",
        annotation=MAGS / "{mag_catalogue}.gff",
    output:
        counts=temp(SUBREAD / "{mag_catalogue}" / "{sample_id}.{library_id}.tsv"),
        summary=SUBREAD / "{mag_catalogue}" / "{sample_id}.{library_id}.tsv.summary",
    log:
        SUBREAD / "{mag_catalogue}" / "{sample_id}.{library_id}.log",
    conda:
        "../../environments/subread.yml"
    params:
        sample_library=lambda w: f"{w.sample_id}.{w.library_id}",
        tmp_out=lambda w: SUBREAD
        / w.mag_catalogue
        / f"{w.sample_id}.{w.library_id}.tsv",
    resources:
        mem_mb=4 * 1024,
    shell:
        """
        featureCounts \
            -F GFF \
            -t CDS \
            -g ID \
            -p \
            -a {input.annotation} \
            -o {params.tmp_out} \
            {input.bam} \
        2> {log} 1>&2
        """


rule quantify__subread__parse_counts:
    input:
        SUBREAD / "{mag_catalogue}" / "{sample_id}.{library_id}.tsv",
    output:
        temp(SUBREAD / "{mag_catalogue}" / "{sample_id}.{library_id}.clean.tsv"),
    log:
        SUBREAD / "{mag_catalogue}" / "{sample_id}.{library_id}.clean.log",
    conda:
        "base"
    shell:
        """
        ( grep -v ^# {input} \
        | cut -f 1,7 \
        > {output} \
        ) 2> {log}
        """


rule quantify__subread__join:
    input:
        lambda w: [
            SUBREAD / w.mag_catalogue / f"{sample_id}.{library_id}.clean.tsv"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        SUBREAD / "{mag_catalogue}.tsv.gz",
    log:
        SUBREAD / "{mag_catalogue}.log",
    params:
        subcommand="join",
    wrapper:
        "v5.2.1/utils/csvtk"


rule quantify__subread__all:
    input:
        [SUBREAD / f"{mag_catalogue}.tsv.gz" for mag_catalogue in MAG_CATALOGUES],
