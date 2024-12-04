rule quantify__multiqc:
    """Collect all reports for the quantify step"""
    input:
        samtools_stats=[
            QUANT_BOWTIE2 / mag_catalogue / f"{sample_id}.{library_id}.stats.tsv"
            for mag_catalogue in MAG_CATALOGUES
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        featurecounts=[
            QUANT_SUBREAD
            / f"{mag_catalogue}"
            / f"{sample_id}.{library_id}.featureCounts.summary"
            for mag_catalogue in MAG_CATALOGUES
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        html=RESULTS / "quantify.html",
        zip=RESULTS / "quantify.zip",
    log:
        RESULTS / "quantify.log",
    params:
        extra="--title quantify --dirs --dirs-depth 1 --fullnames --force",
    resources:
        mem_mb=8 * 1024,
    wrapper:
        "v5.2.1/bio/multiqc"


rule quantify__multiqc__all:
    """Collect all per step reports for the pipeline"""
    input:
        rules.quantify__multiqc.output,
