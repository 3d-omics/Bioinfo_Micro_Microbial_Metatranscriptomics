rule reference_hosts_recompress_genome_one:
    """Link the reference genome to the results directory

    Note: STAR requires the genome decompressed.
    """
    input:
        get_host_genome,
    output:
        HOSTS / "{host_name}.fa",
    log:
        HOSTS / "{host_name}.fa.log"
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


rule reference_hosts_recompress_gtf_one:
    """Link the reference annotation to the results directory

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



rule reference_hosts_recompress_genome:
    input:
        [
            HOSTS / f"{host_name}.fa.gz" for host_name in HOST_NAMES
        ]


rule reference_hosts_recompress_gtf:
    input:
        [
            HOSTS / f"{host_name}.gtf.gz" for host_name in HOST_NAMES
        ]


rule reference_hosts:
    input:
        rules.reference_hosts_recompress_genome.input,
        rules.reference_hosts_recompress_gtf.input,
