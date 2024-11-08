include: "mags_functions.smk"


rule quantify__mags__recompress_fa:
    """Recompress the MAGs into the results directory"""
    input:
        fna=get_mags_fasta,
    output:
        fna=MAGS / "{mag_catalogue}.fa.gz",
    log:
        MAGS / "{mag_catalogue}.fa.log",
    conda:
        "../../environments/htslib.yml"
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


rule quantify__mags__annotation_gtf:
    """Link annotation to the mags"""
    input:
        get_mags_gtf,
    output:
        MAGS / "{mag_catalogue}.gtf",
    log:
        MAGS / "{mag_catalogue}.gtf.log",
    conda:
        "../../environments/htslib.yml"
    cache: True
    shell:
        """
        gzip -dc {input} > {output} 2> {log}
        """


rule quantify__mags__all:
    """Prepare all the MAG catalogues"""
    input:
        [
            MAGS / f"{mag_catalogue}.{extension}"
            for mag_catalogue in MAG_CATALOGUES
            for extension in ["fa.gz", "bed6", "gtf"]
        ],
