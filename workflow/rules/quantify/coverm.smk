include: "coverm_functions.smk"


# CoverM Contig
rule quantify__coverm__genome:
    """calculation of MAG-wise coverage"""
    input:
        bam=BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.bam",
    output:
        tsv=COVERM
        / "genome"
        / "{mag_catalogue}"
        / "{method}"
        / "{sample_id}.{library_id}.tsv.gz",
    log:
        COVERM
        / "genome"
        / "{mag_catalogue}"
        / "{method}"
        / "{sample_id}.{library_id}.log",
    conda:
        "../../environments/coverm.yml"
    params:
        method=get_method,
        min_covered_fraction=get_min_covered_fraction,
        separator=get_separator,
    shell:
        """
        ( coverm genome \
            --bam-files {input.bam} \
            --methods {params.method} \
            --separator "{params.separator}" \
            --min-covered-fraction {params.min_covered_fraction} \
            --output-file /dev/stdout \
        | sed '1 s/Genome/sequence_id/' \
        | cut -f 1 -d " " \
        | gzip \
        > {output.tsv}
        ) 2> {log} 1>&2
        """


rule quantify__coverm__genome__aggregate:
    """Join all the results from coverm, for all assemblies and samples at once, but a single method"""
    input:
        get_tsvs_for_assembly_coverm_genome,
    output:
        tsv=COVERM / "genome.{mag_catalogue}.{method}.tsv.gz",
    log:
        COVERM / "genome.{mag_catalogue}.{method}.log",
    params:
        subcommand="join",
    wrapper:
        "v5.2.1/utils/csvtk"


rule quantify__coverm__genome__all:
    """
    Get the MAG-wise count tables
    """
    input:
        [
            COVERM / f"genome.{mag_catalogue}.{method}.tsv.gz"
            for method in COVERM_CONTIG_METHODS
            for mag_catalogue in MAG_CATALOGUES
        ],


# CoverM contig ----
rule quantify__coverm__contig:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.bam",
    output:
        tsv=COVERM
        / "contig"
        / "{mag_catalogue}"
        / "{method}"
        / "{sample_id}.{library_id}.tsv.gz",
    log:
        COVERM
        / "contig"
        / "{mag_catalogue}"
        / "{method}"
        / "{sample_id}.{library_id}.log",
    conda:
        "../../environments/coverm.yml"
    params:
        method=get_method,
    shell:
        """
        ( coverm contig \
            --bam-files {input.bam} \
            --methods {params.method} \
            --proper-pairs-only \
            --output-file /dev/stdout \
        | sed '1 s/Contig/sequence_id/' \
        | cut -f 1 -d " " \
        | gzip \
        > {output.tsv}
        ) 2> {log} 1>&2
        """


rule quantify__coverm__contig_aggregate:
    """Aggregate coverm contig results"""
    input:
        get_tsvs_for_assembly_coverm_contig,
    output:
        tsv=COVERM / "contig.{mag_catalogue}.{method}.tsv.gz",
    log:
        COVERM / "contig.{mag_catalogue}.{method}.log",
    params:
        subcommand="join",
    wrapper:
        "v5.2.1/utils/csvtk"


rule quantify__coverm__contig__all:
    """
    Get the contig-wise count table
    """
    input:
        [
            COVERM / f"contig.{mag_catalogue}.{method}.tsv.gz"
            for method in COVERM_CONTIG_METHODS
            for mag_catalogue in MAG_CATALOGUES
        ],


rule quantify__coverm__all:
    """
    Get the MAG- and contig-wise count tables
    """
    input:
        rules.quantify__coverm__genome__all.input,
        rules.quantify__coverm__contig__all.input,
