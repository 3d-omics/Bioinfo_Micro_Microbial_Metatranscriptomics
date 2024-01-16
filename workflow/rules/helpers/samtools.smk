rule _helpers__samtools__faidx_fagz:
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


rule _helpers__samtools__faidx_fa:
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


rule _helpers__samtools__index_bam:
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


rule _helpers__samtools__flagstats_cram:
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


rule _helpers__samtools__idxstats_cram:
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
