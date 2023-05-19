rule reference_recompress_genome:
    """Extract the fasta.gz on config.yaml into genome.fa,gz with bgzip"""
    input:
        fa_gz=features["reference"]["dna"],
    output:
        fa_gz=REFERENCE / "genome.fa.gz",
    log:
        REFERENCE / "genome.log",
    conda:
        "../envs/samtools.yml"
    threads: 8
    shell:
        """
        (gzip \
            --decompres \
            --stdout {input.fa_gz} \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
            --stdout \
            /dev/stdin \
        > {output.fa_gz} \
        ) 2> {log}
        """


rule reference:
    """Re-bgzip the reference genome and known variants"""
    input:
        rules.reference_recompress_genome.output,
