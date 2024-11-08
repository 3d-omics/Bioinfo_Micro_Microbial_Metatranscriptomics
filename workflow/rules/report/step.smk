rule report__step__preprocess__:
    """Collect all reports for the preprocess step"""
    input:
        fastqc=rules.preprocess__reads__fastqc__all.input,
        fastp=[
            FASTP / f"{sample_id}.{library_id}_fastp.json"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        kraken2=rules.preprocess__kraken2__all.input,
        ribodetector=rules.preprocess__ribodetector__fastqc__all.input,
        star=rules.preprocess__star__report.input,
    output:
        html=REPORT_STEP / "preprocess.html",
    log:
        REPORT_STEP / "preprocess.log",
    conda:
        "__environment__.yml"
    params:
        dir=REPORT_STEP,
    resources:
        mem_mb=8 * 1024,
    shell:
        """
        multiqc \
            --title preprocess \
            --force \
            --filename preprocess \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report__step__quantify__:
    """Collect all reports for the quantify step"""
    input:
        rules.quantify__bowtie2__report.input,
    output:
        html=REPORT_STEP / "quantify.html",
    log:
        REPORT_STEP / "quantify.log",
    conda:
        "__environment__.yml"
    params:
        dir=REPORT_STEP,
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


rule report__step:
    """Collect all per step reports for the pipeline"""
    input:
        rules.report__step__preprocess__.output,
        rules.report__step__quantify__.output,
