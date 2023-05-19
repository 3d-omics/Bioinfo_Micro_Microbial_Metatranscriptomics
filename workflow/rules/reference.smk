rule reference_set_dna:
    input:
        fa=features["host"]["dna"],
    output:
        fa=REFERENCE / "genome.fa",
    log:
        REFERENCE / "genome.log",
    conda:
        "../envs/empty.yml"
    shell:
        "gzip --decompress --stdout {input.fa} > {output.fa} 2> {log}"


rule reference_set_gtf:
    input:
        gtf=features["host"]["gtf"],
    output:
        gtf=REFERENCE / "annotation.gtf",
    log:
        REFERENCE / "annotation.log",
    conda:
        "../envs/empty.yml"
    shell:
        "gzip --decompress --stdout {input.gtf} > {output.gtf}"


rule reference_join_mags:
    input:
        fna=features["mags"],
    output:
        fna=REFERENCE / "mags.fa.gz",
    log:
        REFERENCE / "mags.log",
    conda:
        "../envs/samtools.yml"
    shell:
        """
        (gzip \
            --decompress \
            --stdout \
            {input.fna} \
        | bgzip \
        >  {output.fna} ) \
        2> {log}
        """


rule reference:
    input:
        rules.reference_set_dna.output.fa,
        rules.reference_set_gtf.output.gtf,
