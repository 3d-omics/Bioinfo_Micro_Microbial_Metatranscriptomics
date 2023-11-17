rule reference_set_dna:
    """Link the reference genome to the results directory"""
    input:
        fa=features["host"]["dna"],
    output:
        fa=REFERENCE / "genome.fa",
    log:
        REFERENCE / "genome.log",
    conda:
        "_env.yml"
    shell:
        "gzip --decompress --stdout {input.fa} > {output.fa} 2> {log}"


rule reference_set_gtf:
    """Link the reference annotation to the results directory"""
    input:
        gtf=features["host"]["gtf"],
    output:
        gtf=REFERENCE / "annotation.gtf",
    log:
        REFERENCE / "annotation.log",
    conda:
        "_env.yml"
    shell:
        "gzip --decompress --stdout {input.gtf} > {output.gtf} 2>{log}"
