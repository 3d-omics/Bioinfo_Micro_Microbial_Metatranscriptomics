rule report_library_one:
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


rule report_library_all:
    """Make a MultiQC report for every library"""
    input:
        [
            REPORT_LIBRARY / f"{sample_id}.{library_id}.html"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],


rule report_library:
    """Make all MultiQC reports per library"""
    input:
        rules.report_library_all.input,

localrules:
    report_library_one,
