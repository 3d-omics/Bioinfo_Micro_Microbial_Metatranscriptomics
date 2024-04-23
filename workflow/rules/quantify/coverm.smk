# CoverM Contig
rule _quantify__coverm__genome:
    """calculation of MAG-wise coverage"""
    input:
        cram=BOWTIE2 / "{mag_catalogue}.{sample_id}.{library_id}.cram",
        crai=BOWTIE2 / "{mag_catalogue}.{sample_id}.{library_id}.cram.crai",
        reference=MAGS / "{mag_catalogue}.fa.gz",
        fai=MAGS / "{mag_catalogue}.fa.gz.fai",
    output:
        tsv=COVERM
        / "genome"
        / "{mag_catalogue}"
        / "{method}"
        / "{sample_id}.{library_id}.tsv",
    log:
        COVERM
        / "genome"
        / "{mag_catalogue}"
        / "{method}"
        / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    params:
        method=get_method,
        min_covered_fraction=get_min_covered_fraction,
        separator=get_separator,
    resources:
        runtime=24 * 60,
        mem_mb=32 * 1024,
    shell:
        """
        ( samtools view \
            --with-header \
            --reference {input.reference} \
            {input.cram} \
        | coverm genome \
            --bam-files /dev/stdin \
            --methods {params.method} \
            --separator "{params.separator}" \
            --min-covered-fraction {params.min_covered_fraction} \
            --output-file {output.tsv} \
        ) 2> {log} 1>&2
        """


rule _quantify__coverm__genome_aggregate:
    """Join all the results from coverm, for all assemblies and samples at once, but a single method"""
    input:
        get_tsvs_for_assembly_coverm_genome,
    output:
        tsv=COVERM / "genome.{mag_catalogue}.{method}.tsv",
    log:
        COVERM / "genome.{mag_catalogue}.{method}.log",
    conda:
        "__environment__.yml"
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


rule quantify__coverm__genome:
    """
    Get the MAG-wise count tables
    """
    input:
        [
            COVERM / f"genome.{mag_catalogue}.{method}.tsv"
            for method in COVERM_CONTIG_METHODS
            for mag_catalogue in MAG_CATALOGUES
        ],


# CoverM contig ----
rule _quantify__coverm__contig:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        cram=BOWTIE2 / "{mag_catalogue}.{sample_id}.{library_id}.cram",
        crai=BOWTIE2 / "{mag_catalogue}.{sample_id}.{library_id}.cram.crai",
        reference=MAGS / "{mag_catalogue}.fa.gz",
        fai=MAGS / "{mag_catalogue}.fa.gz.fai",
    output:
        tsv=COVERM
        / "contig"
        / "{mag_catalogue}"
        / "{method}"
        / "{sample_id}.{library_id}.tsv",
    log:
        COVERM
        / "contig"
        / "{mag_catalogue}"
        / "{method}"
        / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    params:
        method=get_method,
    shell:
        """
        ( samtools view \
            --with-header \
            --reference {input.reference} \
            {input.cram} \
        | coverm contig \
            --bam-files /dev/stdin \
            --methods {params.method} \
            --proper-pairs-only \
            --output-file {output.tsv} \
        ) 2> {log} 1>&2
        """


rule _quantify__coverm__contig_aggregate:
    """Aggregate coverm contig results"""
    input:
        get_tsvs_for_assembly_coverm_contig,
    output:
        tsv=COVERM / "contig.{mag_catalogue}.{method}.tsv",
    log:
        COVERM / "contig.{mag_catalogue}.{method}.log",
    conda:
        "__environment__.yml"
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


rule quantify__coverm__contig:
    """
    Get the contig-wise count table
    """
    input:
        [
            COVERM / f"contig.{mag_catalogue}.{method}.tsv"
            for method in COVERM_CONTIG_METHODS
            for mag_catalogue in MAG_CATALOGUES
        ],


rule quantify__coverm:
    """
    Get the MAG- and contig-wise count tables
    """
    input:
        rules.quantify__coverm__genome.input,
        rules.quantify__coverm__contig.input,
