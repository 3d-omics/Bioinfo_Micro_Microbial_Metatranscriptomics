
rule quantify__multiqc:
    """Collect all reports for the quantify step"""
    input:
        samtools_stats=[
            BOWTIE2 / mag_catalogue / f"{sample_id}.{library_id}.stats.tsv"
            for mag_catalogue in MAG_CATALOGUES
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        html=RESULTS / "quantify.html",
    log:
        RESULTS / "quantify.log",
    conda:
        "../../environments/multiqc.yml"
    params:
        dir=RESULTS,
    resources:
        mem_mb=8 * 1024,
    shell:
        """
        multiqc \
            --title quantify \
            --force \
            --filename quantify \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule quantify__multiqc__all:
    """Collect all per step reports for the pipeline"""
    input:
        rules.quantify__multiqc.output,
