rule quantify__subread__feature_counts:
    input:
        samples=lambda w: [
            QUANT_BOWTIE2 / w.mag_catalogue / f"{sample_id}.{library_id}.bam"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        annotation=QUANT_MAGS / "{mag_catalogue}.gff",
    output:
        multiext(
            str(QUANT_SUBREAD / "{mag_catalogue}"),
            ".featureCounts",
            ".featureCounts.summary",
        ),
    log:
        QUANT_SUBREAD / "{mag_catalogue}.log",
    params:
        extra="-F GFF -t CDS,tRNA,rRNA -g ID -p",
    resources:
        mem_mb=4 * 1024,
    wrapper:
        "v5.2.1/bio/subread/featurecounts"


rule quantify__subread__format:
    input:
        QUANT_SUBREAD / "{mag_catalogue}.featureCounts",
    output:
        QUANT_SUBREAD / "{mag_catalogue}.tsv.gz",
    log:
        QUANT_SUBREAD / "{mag_catalogue}.tsv.log",
    conda:
        "base"
    params:
        workdir=lambda w: str(QUANT_BOWTIE2 / w.mag_catalogue).replace("/", "\\/"),
    shell:
        """
        ( grep -v ^# {input} \
        | cut -f 1,7- \
        | awk '{{ gsub("Geneid", "gene_id", $1); print }}' \
        | awk '{{ gsub("{params.workdir}/", "", $0); print }}' \
        | awk '{{ gsub(".bam", "", $0); print }}' \
        | gzip \
        > {output} \
        ) 2>> {log}
        """


rule quantify__subread__all:
    input:
        [QUANT_SUBREAD / f"{mag_catalogue}.tsv.gz" for mag_catalogue in MAG_CATALOGUES],
