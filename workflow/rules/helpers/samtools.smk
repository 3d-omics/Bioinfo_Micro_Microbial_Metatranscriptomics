rule _fa_gz_fai:
    input:
        "{prefix}.fa.gz",
    output:
        "{prefix}.fa.gz.fai",
    log:
        "{prefix}.fa.fai.gz.log",
    conda:
        "_env.yml"
    shell:
        "samtools faidx {input} 2> {log} 1>&2"


rule _fa_fai:
    input:
        "{prefix}.fa",
    output:
        "{prefix}.fa.fai",
    log:
        "{prefix}.fa.fai.log",
    conda:
        "_env.yml"
    shell:
        "samtools faidx {input} 2> {log} 1>&2"


rule _bai:
    """Generate a bam index"""
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    log:
        "{prefix}.bam.bai.log",
    conda:
        "_env.yml"
    shell:
        "samtools index {input} 2> {log} 1>&2"


rule _crai:
    """Generate a cram index"""
    input:
        "{prefix}.cram",
    output:
        "{prefix}.cram.crai",
    log:
        "{prefix}.cram.crai.log",
    conda:
        "_env.yml"
    shell:
        "samtools index {input} 2> {log} 1>&2"


rule _samtools_flagstats_cram:
    """Compute flagstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        txt="{prefix}.flagstats.txt",
    log:
        "{prefix}.flagstats.log",
    conda:
        "_env.yml"
    shell:
        "samtools flagstats {input.cram} > {output.txt} 2> {log}"


rule _samtools_idxstats_cram:
    """Compute idxstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        tsv="{prefix}.idxstats.tsv",
    log:
        "{prefix}.idxstats.log",
    conda:
        "_env.yml"
    shell:
        "samtools idxstats {input.cram} > {output.tsv} 2> {log}"
