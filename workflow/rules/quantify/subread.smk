rule quantify__subread__feature_counts:
    input:
        bam=BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.bam",
        annotation=MAGS / "{mag_catalogue}.gff",
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
        featureCounts \
            -F GFF \
            -t CDS \
            -g ID \
            -p \
            -a {input.annotation} \
            -o {params.tmp_out} \
            {input.bam} \
        2> {log} 1>&2

        ( grep -v ^# {params.tmp_out} \
        | cut -f 1,7 \
        | gzip \
        > {output.counts} \
        ) 2>> {log}

        rm -rfv {params.tmp_out} 2>> {log} 1>&2
        """


rule quantify__subread__join:
    input:
        lambda w: [
            SUBREAD / w.mag_catalogue / f"{sample_id}.{library_id}.tsv.gz"
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
