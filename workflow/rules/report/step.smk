rule report_step_reads:
    """Collect all reports for the reads step"""
    input:
        rules.reads_fastqc_all.input,
    output:
        html=REPORT_STEP / "reads.html",
    log:
        REPORT_STEP / "reads.log",
    conda:
        "_env.yml"
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


rule report__step__fastp:
    """Collect all reports for the fastp step"""
    input:
        rules.preprocess__fastp__report.input,
    output:
        html=REPORT_STEP / "fastp.html",
    log:
        REPORT_STEP / "fastp.log",
    conda:
        "_env.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title fastp \
            --force \
            --filename fastp \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report__step__ribodetector:
    """Collect all reports for the ribodetector step"""
    input:
        rules.preprocess__ribodetector__fastqc.input,
    output:
        html=REPORT_STEP / "ribodetector.html",
    log:
        REPORT_STEP / "ribodetector.log",
    conda:
        "_env.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title ribodetector \
            --force \
            --filename ribodetector \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report__step__kraken2:
    """Collect all reports for the fastp step"""
    input:
        rules.preprocess__kraken2__report.input,
    output:
        html=REPORT_STEP / "kraken2.html",
    log:
        REPORT_STEP / "kraken2.log",
    conda:
        "_env.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title kraken2 \
            --force \
            --filename kraken2 \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report__step__star:
    """Collect all reports for the star step"""
    input:
        rules.preprocess__star__report.input,
    output:
        html=REPORT_STEP / "star.html",
    log:
        REPORT_STEP / "star.log",
    conda:
        "_env.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title star \
            --force \
            --filename star \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report_step_bowtie2:
    """Collect all reports for the bowtie2 step"""
    input:
        rules.quantification_bowtie2_report_all.input,
    output:
        html=REPORT_STEP / "bowtie2.html",
    log:
        REPORT_STEP / "bowtie2.log",
    conda:
        "_env.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title bowtie2 \
            --force \
            --filename bowtie2 \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report_step:
    """Collect all per step reports for the pipeline"""
    input:
        rules.report_step_reads.output,
        rules.report__step__fastp.output,
        rules.report__step__kraken2.output if features["kraken2_databases"] else [],
        rules.report__step__star.output if features["hosts"] else [],  # No point if no hosts
        rules.report_step_bowtie2.output,
