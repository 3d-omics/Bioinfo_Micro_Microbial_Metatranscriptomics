include: "bowtie2_functions.smk"


rule quantify__bowtie2__build:
    """Build bowtie2 index for the mags"""
    input:
        ref=QUANT_MAGS / "{mag_catalogue}.fa.gz",
    output:
        multiext(
            str(QUANT_INDEX / "{mag_catalogue}"),
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    log:
        QUANT_INDEX / "{mag_catalogue}.log",
    conda:
        "../../environments/bowtie2.yml"
    params:
        extra=params["quantify"]["bowtie2"]["extra"],
    retries: 5
    threads: 24
    resources:
        mem_mb=double_ram(32 * 1024),
        runtime=6 * 60,
    wrapper:
        "v5.2.1/bio/bowtie2/build"


rule quantify__bowtie2__build__all:
    """Build all the bowtie2 indexes"""
    input:
        [
            QUANT_INDEX / f"{mag_catalogue}.{extension}"
            for mag_catalogue in MAG_CATALOGUES
            for extension in [
                "1.bt2",
                "2.bt2",
                "3.bt2",
                "4.bt2",
                "rev.1.bt2",
                "rev.2.bt2",
            ]
        ],


rule quantify__bowtie2__map:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        forward_=PRE_CLEAN / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_CLEAN / "{sample_id}.{library_id}_2.fq.gz",
        bowtie2_index=multiext(
            str(QUANT_INDEX / "{mag_catalogue}"),
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        bam=QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.bam",
    log:
        QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.log",
    conda:
        "../../environments/bowtie2.yml"
    params:
        extra=params["quantify"]["bowtie2"]["extra"],
        samtools_mem=params["quantify"]["bowtie2"]["samtools"]["mem_per_thread"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
        prefix=lambda w: QUANT_INDEX / w.mag_catalogue,
    threads: 24
    resources:
        mem_mb=double_ram(32 * 1024),
        runtime=6 * 60,
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
            -o {output.bam} \
            --threads {threads} \
            --write-index \
        ) 2> {log} 1>&2
        """


rule quantify__bowtie2__map__all:
    """Collect the results of `bowtie2_map_one` for all libraries"""
    input:
        [
            QUANT_BOWTIE2 / f"{mag_catalogue}" / f"{sample_id}.{library_id}.bam"
            for sample_id, library_id in SAMPLE_LIBRARY
            for mag_catalogue in MAG_CATALOGUES
        ],


rule quantify__bowtie2__all:
    """Run bowtie2 on all libraries and generate reports"""
    input:
        rules.quantify__bowtie2__build__all.input,
        rules.quantify__bowtie2__map__all.input,
