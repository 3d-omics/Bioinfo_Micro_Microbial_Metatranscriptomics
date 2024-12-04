include: "fastp_functions.smk"


rule preprocess__fastp:
    """Run fastp on one library"""
    input:
        sample=[
            PRE_READS / "{sample_id}.{library_id}_1.fq.gz",
            PRE_READS / "{sample_id}.{library_id}_2.fq.gz",
        ],
    output:
        trimmed=temp(
            [
                PRE_FASTP / "{sample_id}.{library_id}_1.fq.gz",
                PRE_FASTP / "{sample_id}.{library_id}_2.fq.gz",
            ]
        ),
        html=PRE_FASTP / "{sample_id}.{library_id}.html",
        json=PRE_FASTP / "{sample_id}.{library_id}_fastp.json",
    log:
        PRE_FASTP / "{sample_id}.{library_id}.log",
    params:
        adapters=compose_adapters,
        extra=params["preprocess"]["fastp"]["extra"],
    # group:
    #     "preprocess__{sample_id}.{library_id}"
    threads: 24
    resources:
        mem_mb=4 * 1024,
        runtime=1 * 60,
    wrapper:
        "v5.0.1/bio/fastp"


rule preprocess__fastp__all:
    """Run fastp over all libraries"""
    input:
        [
            PRE_FASTP / f"{sample_id}.{library_id}_{end}.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in [1, 2]
        ],
