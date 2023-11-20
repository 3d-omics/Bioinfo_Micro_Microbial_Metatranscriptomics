rule _reference__mags__recompress:
    """Recpmpress the MAGs into the results directory"""
    input:
        fna=get_mag_catalogue,
    output:
        fna=MAGS / "{mag_catalogue}.fa.gz",
    log:
        MAGS / "{mag_catalogue}.log",
    conda:
        "_env.yml"
    threads: 24
    shell:
        """
        ( gzip \
            --decompress \
            --stdout \
            {input.fna} \
        | bgzip \
            -@ {threads} \
            -l 9 \
        > {output.fna} \
        ) 2> {log}
        """


rule reference__mags:
    input:
        [MAGS / f"{mag_catalogue}.fa.gz" for mag_catalogue in MAG_CATALOGUES],
