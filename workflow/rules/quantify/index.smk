rule quantify__index___:
    """Build bowtie2 index for the mags"""
    input:
        reference=MAGS / "{mag_catalogue}.fa.gz",
        fai=MAGS / "{mag_catalogue}.fa.gz.fai",
    output:
        multiext(
            str(BOWTIE2_INDEX / "{mag_catalogue}"),
            ".1.bt2l",
            ".2.bt2l",
            ".3.bt2l",
            ".4.bt2l",
            ".rev.1.bt2l",
            ".rev.2.bt2l",
        ),
    log:
        BOWTIE2_INDEX / "{mag_catalogue}.log",
    conda:
        "__environment__.yml"
    params:
        extra=params["quantify"]["bowtie2"]["extra"],
        prefix=lambda w: BOWTIE2_INDEX / w.mag_catalogue,
    retries: 5
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            --large-index \
            {params.extra} \
            {input.reference} \
            {params.prefix} \
        2> {log} 1>&2
        """


rule quantify__bowtie2__build:
    """Build all the bowtie2 indexes"""
    input:
        [BOWTIE2_INDEX / f"{mag_catalogue}" for mag_catalogue in MAG_CATALOGUES],
