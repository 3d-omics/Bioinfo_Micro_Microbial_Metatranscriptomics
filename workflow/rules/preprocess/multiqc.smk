rule preprocess__multiqc:
    """Collect all reports for the preprocess step"""
    input:
        fastqc=rules.preprocess__reads__fastqc__all.input,
        fastp=[
            FASTP / f"{sample_id}.{library_id}_fastp.json"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        kraken2=rules.preprocess__kraken2__all.input,
        ribodetector=rules.preprocess__ribodetector__fastqc__all.input,
        star=[
            STAR / host_name / f"{sample_id}.{library_id}.Log.final.out"
            for sample_id, library_id in SAMPLE_LIBRARY
            for host_name in HOST_NAMES
        ],
    output:
        html=RESULTS / "preprocess.html",
    log:
        RESULTS / "preprocess.log",
    conda:
        "../../environments/multiqc.yml"
    params:
        dir=RESULTS,
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
