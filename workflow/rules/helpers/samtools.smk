rule fa_gz_fai:
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


rule fa_fai:
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


rule bai:
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


rule crai:
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


rule dict_fa:
    """Generate a dictionary from a .fa"""
    input:
        "{prefix}.fa",
    output:
        "{prefix}.dict",
    log:
        "{prefix}.dict.log",
    conda:
        "_env.yml"
    shell:
        "samtools dict {input} --output {output} 2> {log} 1>&2"


rule dict_fagz:
    """Generate a dictionary from a fa.gz"""
    input:
        "{prefix}.fa.gz",
    output:
        "{prefix}.dict",
    log:
        "{prefix}.dict.log",
    conda:
        "_env.yml"
    shell:
        "samtools dict {input} --output {output} 2> {log} 1>&2"


# rule samtools_stats_bam:
#     """Compute stats for a bam"""
#     input:
#         bam="{prefix}.bam",
#         bai="{prefix}.bam.bai",
#     output:
#         tsv="{prefix}.stats.tsv",
#     log:
#         "{prefix}.stats.log",
#     conda:
#         "_env.yml"
#     shell:
#         "samtools stats --reference {input.reference} {input.bam} > {output.tsv} 2> {log}"


# rule samtools_flagstats_bam:
#     """Compute flagstats for a bam"""
#     input:
#         bam="{prefix}.bam",
#         bai="{prefix}.bam.bai",
#     output:
#         txt="{prefix}.flagstats.txt",
#     log:
#         "{prefix}.flagstats.log",
#     conda:
#         "_env.yml"
#     shell:
#         "samtools flagstats {input.bam} > {output.txt} 2> {log}"


rule samtools_flagstats_cram:
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


# rule samtools_idxstats_bam:
#     """Compute idxstats for a bam"""
#     input:
#         bam="{prefix}.bam",
#         bai="{prefix}.bam.bai",
#     output:
#         tsv="{prefix}.idxstats.tsv",
#     log:
#         "{prefix}.idxstats.log",
#     conda:
#         "_env.yml"
#     shell:
#         "samtools idxstats {input.bam} > {output.tsv} 2> {log}"


rule samtools_idxstats_cram:
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
