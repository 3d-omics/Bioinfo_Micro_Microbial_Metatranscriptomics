include: "fastp_functions.smk"


rule preprocess__fastp:
    """Run fastp on one library"""
    input:
        sample=[
            READS / "{sample_id}.{library_id}_1.fq.gz",
            READS / "{sample_id}.{library_id}_2.fq.gz",
        ],
    output:
        trimmed=[
            FASTP / "{sample_id}.{library_id}_1.fq.gz",
            FASTP / "{sample_id}.{library_id}_2.fq.gz",
        ],
        html=FASTP / "{sample_id}.{library_id}.html",
        json=FASTP / "{sample_id}.{library_id}_fastp.json",
    log:
        FASTP / "{sample_id}.{library_id}.log",
    params:
        adapters=compose_adapters,
        extra=params["preprocess"]["fastp"]["extra"],
    wrapper:
        "v5.0.1/bio/fastp"


rule preprocess__fastp__all:
    """Run fastp over all libraries"""
    input:
        [
            FASTP / f"{sample_id}.{library_id}_{end}.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in [1, 2]
        ],
