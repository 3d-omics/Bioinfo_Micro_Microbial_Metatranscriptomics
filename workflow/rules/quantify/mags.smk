rule quantify__mags__recompress_fa:
    """Recompress the MAGs into the results directory"""
    input:
        lambda w: features["mag_catalogues"][w.mag_catalogue]["fasta"],
    output:
        MAGS / "{mag_catalogue}.fa.gz",
    log:
        MAGS / "{mag_catalogue}.fa.log",
    conda:
        "../../environments/htslib.yml"
    cache: True
    threads: 8
    shell:
        """
        ( gzip \
            --decompress \
            --stdout \
            {input} \
        | bgzip \
            -@ {threads} \
        > {output} \
        ) 2> {log}
        """


rule quantify__mags__annotation_gff:
    """Link annotation to the mags"""
    input:
        lambda w: features["mag_catalogues"][w.mag_catalogue]["gff"],
    output:
        MAGS / "{mag_catalogue}.gff",
    log:
        MAGS / "{mag_catalogue}.gff.log",
    conda:
        "base"
    cache: True
    shell:
        """
        gzip --decompress --stdout {input} > {output} 2> {log}
        """


rule quantify__mags__all:
    """Prepare all the MAG catalogues"""
    input:
        [
            MAGS / f"{mag_catalogue}.{extension}"
            for mag_catalogue in MAG_CATALOGUES
            for extension in ["fa.gz", "gff"]
        ],
