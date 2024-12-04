rule quantify__subread__feature_counts:
    input:
        samples=QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.bam",
        annotation=QUANT_MAGS / "{mag_catalogue}.gff",
    output:
        multiext(
            str(QUANT_SUBREAD / "{mag_catalogue}" / "{sample_id}.{library_id}"),
            ".featureCounts",
            ".featureCounts.summary",
        ),
    log:
        QUANT_SUBREAD / "{mag_catalogue}" / "{sample_id}.{library_id}.log",
    # conda:
    #     "../../environments/subread.yml"
    params:
        extra="-F GFF -t CDS,tRNA,rRNA -g ID -p",
    resources:
        mem_mb=4 * 1024,
    wrapper:
        "v5.2.1/bio/subread/featurecounts"


rule quantify__subread__extract:
    input:
        QUANT_SUBREAD / "{mag_catalogue}" / "{sample_id}.{library_id}.featureCounts",
    output:
        QUANT_SUBREAD / "{mag_catalogue}" / "{sample_id}.{library_id}.tsv",
    log:
        QUANT_SUBREAD / "{mag_catalogue}" / "{sample_id}.{library_id}.tsv.log",
    conda:
        "base"
    params:
        sample_library=lambda w: f"{w.sample_id}_{w.library_id}",
    shell:
        """
        ( grep -v ^# {input} \
        | cut -f 1,7 \
        | awk '$0 = NR == 1 ? replace : $0' replace='gene_id\t{params.sample_library}' \
        > {output} \
        ) 2>> {log}
        """


rule quantify__subread__join:
    input:
        lambda w: [
            QUANT_SUBREAD / w.mag_catalogue / f"{sample_id}.{library_id}.tsv"
            for sample_id, library_id in SAMPLE_LIBRARY
        ]
        + ["/dev/null"],
    output:
        QUANT_SUBREAD / "{mag_catalogue}.tsv.gz",
    log:
        QUANT_SUBREAD / "{mag_catalogue}.log",
    params:
        subcommand="join",
        extra="--left-join --tabs --out-tabs",
        col1="gene_id",
        col2="gene_id",
    resources:
        mem_mb=double_ram(32 * 1024),
    retries: 5
    wrapper:
        "v5.2.1/utils/csvtk"


rule quantify__subread__all:
    input:
        [QUANT_SUBREAD / f"{mag_catalogue}.tsv.gz" for mag_catalogue in MAG_CATALOGUES],
