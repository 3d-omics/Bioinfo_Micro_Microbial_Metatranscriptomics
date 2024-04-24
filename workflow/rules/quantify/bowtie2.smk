rule quantify__bowtie2__map__:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        forward_=get_forward_for_bowtie2,
        reverse_=get_reverse_for_bowtie2,
        bowtie2_index=multiext(
            str(BOWTIE2_INDEX / "{mag_catalogue}"),
            ".1.bt2l",
            ".2.bt2l",
            ".3.bt2l",
            ".4.bt2l",
            ".rev.1.bt2l",
            ".rev.2.bt2l",
        ),
        reference=MAGS / "{mag_catalogue}.fa.gz",
        fai=MAGS / "{mag_catalogue}.fa.gz.fai",
    output:
        cram=BOWTIE2 / "{mag_catalogue}.{sample_id}.{library_id}.cram",
        crai=BOWTIE2 / "{mag_catalogue}.{sample_id}.{library_id}.cram.crai",
    log:
        BOWTIE2 / "{mag_catalogue}.{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    params:
        extra=params["quantify"]["bowtie2"]["extra"],
        samtools_mem=params["quantify"]["bowtie2"]["samtools"]["mem_per_thread"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
        prefix=lambda w: BOWTIE2_INDEX / w.mag_catalogue,
    shell:
        """
        ( bowtie2 \
            -x {params.prefix} \
            -1 {input.forward_} \
            -2 {input.reverse_} \
            --threads {threads} \
            --rg-id '{params.rg_id}' \
            --rg '{params.rg_extra}' \
            {params.extra} \
        | samtools sort \
            -m {params.samtools_mem} \
            -o {output.cram} \
            --reference {input.reference} \
            --threads {threads} \
            --write-index \
        ) 2> {log} 1>&2
        """


rule quantify__bowtie2__map:
    """Collect the results of `bowtie2_map_one` for all libraries"""
    input:
        [
            BOWTIE2 / f"{mag_catalogue}.{sample_id}.{library_id}.cram"
            for sample_id, library_id in SAMPLE_LIBRARY
            for mag_catalogue in MAG_CATALOGUES
        ],


rule quantify__bowtie2__report:
    """Generate bowtie2 reports for all libraries:
    - samtools stats
    - samtools flagstats
    - samtools idxstats
    """
    input:
        [
            BOWTIE2 / f"{mag_catalogue}.{sample_id}.{library_id}.{report}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for report in BAM_REPORTS
            for mag_catalogue in MAG_CATALOGUES
        ],


rule quantify__bowtie2:
    """Run bowtie2 on all libraries and generate reports"""
    input:
        rules.quantify__bowtie2__map.input,
        rules.quantify__bowtie2__report.input,
