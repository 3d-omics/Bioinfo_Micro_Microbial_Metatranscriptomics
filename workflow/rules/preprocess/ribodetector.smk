rule preprocess__ribodetector__cpu:
    """Run ribodetector on one library

    ribodetector filters out rRNA reads from a library
    """
    input:
        forward_=FASTP / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=FASTP / "{sample_id}.{library_id}_2.fq.gz",
    output:
        forward_=RIBODETECTOR / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=RIBODETECTOR / "{sample_id}.{library_id}_2.fq.gz",
    log:
        RIBODETECTOR / "{sample_id}.{library_id}.log",
    conda:
        "../../environments/ribodetector.yml"
    params:
        average_length=params["preprocess"]["ribodetector"]["average_length"],
        chunk_size=params["preprocess"]["ribodetector"]["chunk_size"],
    shell:
        """
        ribodetector_cpu \
            --input \
                {input.forward_} \
                {input.reverse_} \
            --output \
                {output.forward_} \
                {output.reverse_} \
            --len {params.average_length} \
            --ensure rrna \
            --threads {threads} \
            --chunk_size {params.chunk_size} \
        2> {log} 1>&2
        """


rule preprocess__ribodetector__cpu__all:
    """Run ribodetector_find_one over all libraries"""
    input:
        [
            RIBODETECTOR / f"{sample_id}.{library_id}_{end}.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in [1, 2]
        ],


rule preprocess__ribodetector__fastqc__all:
    input:
        [
            RIBODETECTOR / f"{sample_id}.{library_id}_{end}_fastqc.zip"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in [1, 2]
        ],


rule preprocess__ribodetector__all:
    input:
        rules.preprocess__ribodetector__cpu__all.input,
        rules.preprocess__ribodetector__fastqc__all.input,
