rule preprocess__hosts__decompress_fa:
    """Decompress the genome to its place

    Note: STAR requires the genome decompressed.
    """
    input:
        lambda w: features["hosts"][w.host_name]["genome"],
    output:
        PRE_HOSTS / "{host_name}.fa",
    log:
        PRE_HOSTS / "{host_name}.fa.log",
    cache: True
    conda:
        "base"
    shell:
        """
        gzip \
            --decompress \
            --stdout \
            {input} \
        > {output} \
        2> {log}
        """


rule preprocess__hosts__decompress_gtf:
    """Decomplress the GTF annotation to its place

    Note: STAR requires the annotation decompressed
    """
    input:
        lambda w: features["hosts"][w.host_name]["gtf"],
    output:
        PRE_HOSTS / "{host_name}.gtf",
    log:
        PRE_HOSTS / "{host_name}.gtf.log",
    cache: True
    conda:
        "base"
    shell:
        """
        gzip \
            --decompress \
            --stdout \
            {input} \
        > {output} \
        2>{log}
        """


rule preprocess__hosts__decompress_fa__all:
    """Decompress all the hosts fasta files"""
    input:
        [PRE_HOSTS / f"{host_name}.fa" for host_name in HOST_NAMES],


rule preprocess__hosts__decompress_gtf__all:
    """Decompress all the host GTF annotations"""
    input:
        [PRE_HOSTS / f"{host_name}.gtf" for host_name in HOST_NAMES],


rule preprocess__hosts__all:
    """Prepare all the host files"""
    input:
        rules.preprocess__hosts__decompress_fa__all.input,
        rules.preprocess__hosts__decompress_gtf__all.input,
