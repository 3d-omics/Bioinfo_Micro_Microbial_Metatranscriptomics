include: "clean_functions.smk"


rule preprocess__clean:
    input:
        forward_=get_final_forward,
        reverse_=get_final_reverse,
    output:
        forward_=CLEAN / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=CLEAN / "{sample_id}.{library_id}_2.fq.gz",
    log:
        CLEAN / "{sample_id}.{library_id}.log",
    conda:
        "../../environments/star.yml"
    group:
        "{sample_id}.{library_id}"
    threads: 24
    resources:
        runtime=1 * 60,
    shell:
        """
        ( gzip \
            --decompress \
            --stdout \
            {input.forward_} \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
        > {output.forward_} \
        ) 2> {log}

        ( gzip \
            --decompress \
            --stdout \
            {input.reverse_} \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
        > {output.reverse_} \
        ) 2>> {log}
        """


rule preprocess__clean__all:
    input:
        [
            CLEAN / f"{sample_id}.{library_id}_{end}.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in [1, 2]
        ],
