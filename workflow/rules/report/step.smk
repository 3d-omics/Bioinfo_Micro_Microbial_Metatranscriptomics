rule report__step__reads:
    """Collect all reports for the reads step"""
    input:
        rules.reads__fastqc.input,
    output:
        html=REPORT_STEP / "reads.html",
    log:
        REPORT_STEP / "reads.log",
    conda:
        "__environment__.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --filename reads \
            --title reads \
            --force \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report__step__preprocess:
    """Collect all reports for the preprocess step"""
    input:
        rules.preprocess__fastp__report.input,
        rules.preprocess__kraken2__report.input,
        rules.preprocess__ribodetector__fastqc.input,
        rules.preprocess__star__report.input,
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


rule report__step__quantify:
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
        rules.report__step__reads.output,
        rules.report__step__preprocess.output,
        # rules.report__step__kraken2.output if features["kraken2_databases"] else [],
        # rules.report__step__star.output if features["hosts"] else [],  # No point if no hosts
        rules.report__step__quantify.output,
