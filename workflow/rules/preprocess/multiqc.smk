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
        zip=RESULTS / "preprocess.zip",
    log:
        RESULTS / "preprocess.log",
    params:
        extra="--title preprocess --dirs --dirs-depth 1 --fullnames --force",
    resources:
        mem_mb=8 * 1024,
    wrapper:
        "v5.2.1/bio/multiqc"
