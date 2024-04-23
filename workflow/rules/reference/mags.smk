rule reference__mags__:
    """Recompress the MAGs into the results directory"""
    input:
        fna=get_mag_catalogue,
    output:
        fna=MAGS / "{mag_catalogue}.fa.gz",
    log:
        MAGS / "{mag_catalogue}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        ( gzip \
            --decompress \
            --stdout \
            {input.fna} \
        | bgzip \
            -@ {threads} \
        > {output.fna} \
        ) 2> {log}
        """


rule reference__mags:
    """Prepare all the MAG catalogues"""
    input:
        [MAGS / f"{mag_catalogue}.fa.gz" for mag_catalogue in MAG_CATALOGUES],
