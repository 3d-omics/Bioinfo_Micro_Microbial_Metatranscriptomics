rule _reference__hosts__decompress_fa:
    """Decompress the genome to its place

    Note: STAR requires the genome decompressed.
    """
    input:
        get_host_genome,
    output:
        HOSTS / "{host_name}.fa",
    log:
        HOSTS / "{host_name}.fa.log",
    conda:
        "_env.yml"
    shell:
        """
        gzip \
            --decompress \
            --stdout \
            {input} \
        > {output} \
        2> {log}
        """


rule _reference__hosts__decompress_gtf:
    """Decomplress the GTF annotation to its place

    Note: STAR requires the annotation decompressed
    """
    input:
        get_host_annotation,
    output:
        HOSTS / "{host_name}.gtf",
    log:
        HOSTS / "{host_name}.gtf.log",
    conda:
        "_env.yml"
    shell:
        """
        gzip \
            --decompress \
            --stdout \
            {input} \
        > {output} \
        2>{log}
        """


rule reference__hosts__recompress_fa:
    input:
        [HOSTS / f"{host_name}.fa.gz" for host_name in HOST_NAMES],


rule reference__hosts__recompress_gtf:
    input:
        [HOSTS / f"{host_name}.gtf.gz" for host_name in HOST_NAMES],


rule reference__hosts:
    input:
        rules.reference__hosts__recompress_fa.input,
        rules.reference__hosts__recompress_gtf.input,
