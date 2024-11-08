include: "fastp_functions.smk"


rule preprocess__fastp:
    """Run fastp on one library"""
    input:
        forward_=READS / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=READS / "{sample_id}.{library_id}_2.fq.gz",
    output:
        forward_=FASTP / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=FASTP / "{sample_id}.{library_id}_2.fq.gz",
        unpaired1=FASTP / "{sample_id}.{library_id}_u1.fq.gz",
        unpaired2=FASTP / "{sample_id}.{library_id}_u2.fq.gz",
        html=FASTP / "{sample_id}.{library_id}.html",
        json=FASTP / "{sample_id}.{library_id}_fastp.json",
    log:
        FASTP / "{sample_id}.{library_id}.log",
    params:
        adapter_forward=get_forward_adapter,
        adapter_reverse=get_reverse_adapter,
        extra=params["preprocess"]["fastp"]["extra"],
    conda:
        "__environment__.yml"
    shell:
        """
        fastp \
            --in1 {input.forward_} \
            --in2 {input.reverse_} \
            --out1 {output.forward_} \
            --out2 {output.reverse_} \
            --unpaired1 {output.unpaired1} \
            --unpaired2 {output.unpaired2} \
            --html {output.html} \
            --json {output.json} \
            --verbose \
            --adapter_sequence {params.adapter_forward} \
            --adapter_sequence_r2 {params.adapter_reverse} \
            --thread {threads} \
            {params.extra} \
        2> {log}
        """


rule preprocess__fastp__all:
    """Run fastp over all libraries"""
    input:
        [
            FASTP / f"{sample_id}.{library_id}_{end}.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in "1 2 u1 u2".split(" ")
        ],
