rule helpers__samtools__faidx_fagz__:
    """Index a .fa.gz file"""
    input:
        "{prefix}.fa.gz",
    output:
        "{prefix}.fa.gz.fai",
    log:
        "{prefix}.fa.fai.gz.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools faidx {input} 2> {log} 1>&2"


rule helpers__samtools__faidx_fa__:
    """Index a .fa file"""
    input:
        "{prefix}.fa",
    output:
        "{prefix}.fa.fai",
    log:
        "{prefix}.fa.fai.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools faidx {input} 2> {log} 1>&2"


rule helpers__samtools__bai__:
    """Generate a bam index"""
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    log:
        "{prefix}.bam.bai.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools index {input} 2> {log} 1>&2"


rule helpers__samtools__flagstats_cram__:
    """Compute flagstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        txt="{prefix}.flagstats.txt",
    log:
        "{prefix}.flagstats.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools flagstats {input.cram} > {output.txt} 2> {log}"


rule helpers__samtools__idxstats_cram__:
    """Compute idxstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        tsv="{prefix}.idxstats.tsv",
    log:
        "{prefix}.idxstats.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools idxstats {input.cram} > {output.tsv} 2> {log}"
