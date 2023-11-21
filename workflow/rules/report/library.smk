rule _report__library:
    """Make a MultiQC report for a single library"""
    input:
        READS / "{sample_id}.{library_id}_1_fastqc.zip",
        READS / "{sample_id}.{library_id}_2_fastqc.zip",
        FASTP / "{sample_id}.{library_id}_fastp.json",
        FASTP / "{sample_id}.{library_id}_1_fastqc.zip",
        FASTP / "{sample_id}.{library_id}_2_fastqc.zip",
        RIBODETECTOR / "{sample_id}.{library_id}_1_fastqc.zip",
        RIBODETECTOR / "{sample_id}.{library_id}_2_fastqc.zip",
        get_kraken2_for_library_report,
        get_samtools_for_library_report,
        get_star_for_library_report,
    output:
        REPORT_LIBRARY / "{sample_id}.{library_id}.html",
    log:
        REPORT_LIBRARY / "{sample_id}.{library_id}.log",
    conda:
        "_env.yml"
    params:
        library="{sample_id}.{library_id}",
        out_dir=REPORT_LIBRARY,
    shell:
        """
        multiqc \
            --title {params.library} \
            --force \
            --filename {params.library} \
            --outdir {params.out_dir} \
            --dirs \
            --dirs-depth 1 \
            {input} \
        2> {log} 1>&2
        """


rule report__library:
    """Make a MultiQC report for every library"""
    input:
        [
            REPORT_LIBRARY / f"{sample_id}.{library_id}.html"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
