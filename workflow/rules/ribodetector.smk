rule ribodetector_filter_one:
    input:
        forward_=FASTP / "{sample}.{library}_1.fq.gz",
        reverse_=FASTP / "{sample}.{library}_2.fq.gz",
    output:
        forward_=temp(RIBODETECTOR / "{sample}.{library}_1.fq.gz"),
        reverse_=temp(RIBODETECTOR / "{sample}.{library}_2.fq.gz"),
    log:
        RIBODETECTOR / "interleaved/{sample}.{library}.log",
    threads: 24
    params:
        len=100
    conda:
        "../envs/ribodetector.yml",
    shell:
        """
        ribodetector_cpu \
            --len {params.len} \
            --input \
                {input.forward_} \
                {input.reverse_} \
            --output \
                >(pigz --fast > {output.forward_}) \
                >(pigz --fast > {output.reverse_}) \
            --ensure rrna \
            --threads {threads} \
        2> {log} 1>&2
        """


rule ribodetector_find_all:
    input:
        [
            RIBODETECTOR / f"{sample}.{library}_{end}.fq.gz"
            for sample, library  in SAMPLE_LIB
            for end in ["1", "2"]
        ]


rule ribodetector_fastqc_all:
    """Run fastqc over all libraries"""
    input:
        [
            RIBODETECTOR / f"{sample}.{library}_{end}_fastqc.{extension}"
            for sample, library in SAMPLE_LIB
            for end in ["1", "2"]
            for extension in ["html", "zip"]
        ],


rule ribodetector:
    input:
        rules.ribodetector_find_all.input,
        rules.ribodetector_fastqc_all.input,
