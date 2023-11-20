rule quantification_bowtie2_build_one:
    """Build bowtie2 index for the mags"""
    input:
        reference=MAGS / "{mag_catalogue}.fa.gz",
        fai = MAGS / "{mag_catalogue}.fa.gz.fai",
    output:
        prefix = touch(BOWTIE2_INDEX / "{mag_catalogue}")
    log:
        BOWTIE2_INDEX / "{mag_catalogue}.log",
    conda:
        "_env.yml"
    params:
        extra=params["quantification"]["bowtie2"]["extra"],
    threads: 24
    resources:
        mem_mb=double_ram(32),
        runtime=6 * 60,
    retries: 5
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {params.extra} \
            {input.reference} \
            {output.prefix} \
        2> {log} 1>&2
        """


rule quantification_bowtie2_map_one:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        forward_=get_forward_for_bowtie2,
        reverse_=get_reverse_for_bowtie2,
        bowtie2_index=BOWTIE2_INDEX / "{mag_catalogue}",
        reference=MAGS / "{mag_catalogue}.fa.gz",
        fai = MAGS / "{mag_catalogue}.fa.gz.fai",
    output:
        cram=BOWTIE2 / "{mag_catalogue}.{sample_id}.{library_id}.cram",
    log:
        BOWTIE2 / "{mag_catalogue}.{sample_id}.{library_id}.log",
    params:
        extra=params["quantification"]["bowtie2"]["extra"],
        samtools_mem=params["quantification"]["bowtie2"]["samtools"]["mem_per_thread"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    threads: 24
    conda:
        "_env.yml"
    resources:
        mem_mb=double_ram(32),
        runtime=1440,
    shell:
        """
        ( bowtie2 \
            -x {input.bowtie2_index} \
            -1 {input.forward_} \
            -2 {input.reverse_} \
            --threads {threads} \
            --rg-id '{params.rg_id}' \
            --rg '{params.rg_extra}' \
            {params.extra} \
        | samtools sort \
            -l 9 \
            -M \
            -m {params.samtools_mem} \
            -o {output.cram} \
            --reference {input.reference} \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


rule quantification_bowtie2_map_all:
    """Collect the results of `bowtie2_map_one` for all libraries"""
    input:
        [
            BOWTIE2 / f"{mag_catalogue}.{sample_id}.{library_id}.cram"
            for sample_id, library_id in SAMPLE_LIBRARY
            for mag_catalogue in MAG_CATALOGUES
        ],


rule quantification_bowtie2_report_all:
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


rule quantification_bowtie2:
    """Run bowtie2 on all libraries and generate reports"""
    input:
        rules.quantification_bowtie2_map_all.input,
        rules.quantification_bowtie2_report_all.input,
