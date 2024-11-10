rule preprocess__multiqc__all:
    """Collect all reports for the preprocess step"""
    input:
        fastqc=rules.preprocess__reads__fastqc__all.input,
        fastp=[
            FASTP / f"{sample_id}.{library_id}_fastp.json"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        ribodetector=[
            RIBODETECTOR / f"{sample_id}.{library_id}_{end}_fastqc.zip"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in [1, 2]
        ],
        star=[
            STAR / f"{host_name}.{sample_id}.{library_id}.Log.final.out"
            for sample_id, library_id in SAMPLE_LIBRARY
            for host_name in HOST_NAMES
        ],
        clean=[
            CLEAN / f"{sample_id}.{library_id}_{end}_fastqc.zip"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in [1, 2]
        ],
        kraken2=[
            KRAKEN2 / kraken_db / f"{sample_id}.{library_id}.report"
            for sample_id, library_id in SAMPLE_LIBRARY
            for kraken_db in features["databases"]["kraken2"]
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
