# CRAM to BAM
rule coverm_cram_to_bam_one:
    """Convert cram to bam

    Note: this step is needed because coverm probably does not support cram. The
    log from coverm shows failures to get the reference online, but nonetheless
    it works.
    """
    input:
        cram=BOWTIE2 / "{mag_catalogue}.{sample_id}.{library_id}.cram",
        crai=BOWTIE2 / "{mag_catalogue}.{sample_id}.{library_id}.cram.crai",
        reference=MAGS / "{mag_catalogue}.fa.gz",
        fai=MAGS / "{mag_catalogue}.fa.gz.fai",
    output:
        bam=temp(COVERM / "bams" / "{mag_catalogue}.{sample_id}.{library_id}.bam"),
    log:
        COVERM / "bams" / "{mag_catalogue}.{sample_id}.{library_id}.bam.log",
    conda:
        "_env.yml"
    resources:
        runtime=1 * 60,
        mem_mb=4 * 1024,
    shell:
        """
        samtools view \
            -F 4 \
            --reference {input.reference} \
            --output {output.bam} \
            --fast \
            {input.cram} \
        2> {log} 1>&2
        """


rule coverm_cram_to_bam:
    input:
        [
            COVERM / "bams" / f"{mag_catalogue}.{sample_id}.{library_id}.bam"
            for mag_catalogue in MAG_CATALOGUES
            for sample_id, library_id in SAMPLE_LIB
        ]


# CoverM Contig
rule coverm_genome_one:
    """calculation of mag-wise coverage"""
    input:
        bam = COVERM / "bams" / "{mag_catalogue}.{sample_id}.{library_id}.bam",
        bai = COVERM / "bams" / "{mag_catalogue}.{sample_id}.{library_id}.bam.bai",
    output:
        tsv=COVERM / "genome" / "{mag_catalogue}" / "{method}" / "{sample_id}.{library_id}.tsv",
    log:
        COVERM / "genome" / "{mag_catalogue}" / "{method}" / "{sample_id}.{library_id}.log"
    conda:
        "_env.yml"
    params:
        method=get_method,
        min_covered_fraction=params["coverm"]["genome"]["min_covered_fraction"],
        separator=params["coverm"]["separator"],
    resources:
        runtime=24 * 60,
        mem_mb=32 * 1024
    shell:
        """
        coverm genome \
            --bam-files {input.bam} \
            --methods {params.method} \
            --separator "{params.separator}" \
            --min-covered-fraction {params.min_covered_fraction} \
            --output-file {output.tsv} \
        2> {log} 1>&2
        """


rule coverm_aggregate_genome:
    """Join all the results from coverm, for all assemblies and samples at once, but a single method"""
    input:
        get_tsvs_for_assembly_coverm_genome,
    output:
        tsv=COVERM / "genome.{mag_catalogue}.{method}.tsv",
    log:
        COVERM / "genome.{mag_catalogue}.{method}.log",
    conda:
        "_env.yml"
    params:
        input_dir=compose_input_dir_for_coverm_genome_aggregate,
    resources:
        mem_mb=8 * 1024,
    shell:
        """
        Rscript --vanilla workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


# CoverM contig ----
rule coverm_contig_one:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=COVERM / "bams" / "{mag_catalogue}.{sample_id}.{library_id}.bam",
        bai=COVERM / "bams" / "{mag_catalogue}.{sample_id}.{library_id}.bam.bai",
        reference=MAGS / "{mag_catalogue}.fa.gz",
        fai=MAGS / "{mag_catalogue}.fa.gz.fai",
    output:
        tsv=COVERM / "contig" / "{mag_catalogue}" / "{method}" / "{sample_id}.{library_id}.tsv",
    log:
        COVERM / "contig" / "{mag_catalogue}" / "{method}" / "{sample_id}.{library_id}.log",
    conda:
        "_env.yml"
    params:
        method=get_method,
    shell:
        """
        coverm contig \
            --bam-files {input.bam} \
            --methods {params.method} \
            --proper-pairs-only \
            --output-file {output.tsv} \
        2> {log} 1>&2
        """


rule coverm_contig_aggregate:
    """Aggregate coverm contig results"""
    input:
        get_tsvs_for_assembly_coverm_contig,
    output:
        tsv=COVERM / "contig.{mag_catalogue}.{method}.tsv",
    log:
        COVERM / "contig.{mag_catalogue}.{method}.log",
    conda:
        "_env.yml"
    params:
        input_dir=compose_input_dir_for_coverm_contig_aggregate,
    resources:
        mem_mb=8 * 1024,
    shell:
        """
        Rscript --vanilla workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule coverm_contig:
    input:
        [
            COVERM / f"contig.{mag_catalogue}.{method}.tsv"
            for method in COVERM_CONTIG_METHODS
            for mag_catalogue in MAG_CATALOGUES
        ]


rule coverm_genome:
    input:
        [
            COVERM / f"genome.{mag_catalogue}.{method}.tsv"
            for method in COVERM_CONTIG_METHODS
            for mag_catalogue in MAG_CATALOGUES
        ]


rule coverm:
    input:
        rules.coverm_genome.input,
        rules.coverm_contig.input,
