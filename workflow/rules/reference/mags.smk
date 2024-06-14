rule reference__mags__fasta__:
    """Recompress the MAGs into the results directory"""
    input:
        fna=get_mags_fasta,
    output:
        fna=MAGS / "{mag_catalogue}.fa.gz",
    log:
        MAGS / "{mag_catalogue}.fa.log",
    conda:
        "__environment__.yml"
    cache: True
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


rule reference__mags__annotation_bed6__:
    """Link annotation to the mags"""
    input:
        get_mags_bed6,
    output:
        MAGS / "{mag_catalogue}.bed6",
    log:
        MAGS / "{mag_catalogue}.bed6.log",
    conda:
        "__environment__.yml"
    cache: True
    shell:
        """
        bedtools sort -i {input} > {output} 2> {log}
        """


rule reference__mags__annotation_gtf__:
    """Link annotation to the mags"""
    input:
        get_mags_gtf,
    output:
        MAGS / "{mag_catalogue}.gtf",
    log:
        MAGS / "{mag_catalogue}.gtf.log",
    conda:
        "__environment__.yml"
    cache: True
    shell:
        """
        ln --symbolic $(readlink --canonicalize {input}) {output} 2> {log}
        """


rule reference__mags:
    """Prepare all the MAG catalogues"""
    input:
        [MAGS / f"{mag_catalogue}.fa.gz" for mag_catalogue in MAG_CATALOGUES],
        [MAGS / f"{mag_catalogue}.bed6" for mag_catalogue in MAG_CATALOGUES],
        [MAGS / f"{mag_catalogue}.gtf" for mag_catalogue in MAG_CATALOGUES],
