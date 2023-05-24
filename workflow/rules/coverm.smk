rule coverm_overall:
    input:
        crams=[BOWTIE2 / f"{sample}.{library}.cram" for sample, library in SAMPLE_LIB],
        mags=REFERENCE / "mags.fa.gz",
    output:
        COVERM / "coverm_overall.tsv",
    log:
        COVERM / "coverm_overall.log",
    conda:
        "../envs/coverm.yml"
    params:
        methods=params["coverm"]["genome"]["methods"],
        min_covered_fraction=params["coverm"]["genome"]["min_covered_fraction"],
    threads: 24
    shell:
        """
        coverm genome \
            --bam-files {input.crams} \
            --methods {params.methods} \
            --separator _ \
            --threads {threads} \
            --min-covered-fraction {params.min_covered_fraction} \
        > {output} \
        2> {log}
        """


rule coverm_contig:
    input:
        crams=[BOWTIE2 / f"{sample}.{library}.cram" for sample, library in SAMPLE_LIB],
        mags=REFERENCE / "mags.fa.gz",
    output:
        COVERM / "coverm_contig.tsv",
    log:
        COVERM / "coverm_contig.log",
    conda:
        "../envs/coverm.yml"
    params:
        methods=params["coverm"]["contig"]["methods"],
    threads: 24
    shell:
        """
        coverm contig \
            --bam-files {input.crams} \
            --methods {params.methods} \
            --proper-pairs-only \
        > {output} \
        2> {log}
        """


rule coverm:
    input:
        rules.coverm_overall.output,
        rules.coverm_contig.output,
