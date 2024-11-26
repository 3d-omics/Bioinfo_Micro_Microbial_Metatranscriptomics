# CoverM Genome ----
rule quantify__coverm__genome:
    """calculation of MAG-wise coverage"""
    input:
        BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.bam",
    output:
        temp(
            COVERM / "genome" / "{mag_catalogue}.{method}.{sample_id}.{library_id}.tsv"
        ),
    log:
        COVERM / "genome" / "{mag_catalogue}.{method}.{sample_id}.{library_id}.log",
    conda:
        "../../environments/coverm.yml"
    params:
        method=lambda w: w.method,
        min_covered_fraction=params["quantify"]["coverm"]["genome"][
            "min_covered_fraction"
        ],
        separator=params["quantify"]["coverm"]["separator"],
    shell:
        """
        ( coverm genome \
            --bam-files {input} \
            --methods {params.method} \
            --separator "{params.separator}" \
            --min-covered-fraction {params.min_covered_fraction} \
            --output-file /dev/stdout \
        | sed '1 s/Genome/sequence_id/' \
        | cut -f 1 -d " " \
        > {output} \
        ) 2> {log} 1>&2
        """


rule quantify__coverm__genome__join:
    """Join all the results from coverm, for all assemblies and samples at once, but a single method"""
    input:
        lambda w: [
            COVERM
            / "genome"
            / f"{w.mag_catalogue}.{w.method}.{sample_id}.{library_id}.tsv"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        COVERM / "genome.{mag_catalogue}.{method}.tsv.gz",
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
        BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.bam",
    output:
        temp(
            COVERM / "contig" / "{mag_catalogue}.{method}.{sample_id}.{library_id}.tsv"
        ),
    log:
        COVERM / "contig" / "{mag_catalogue}.{method}.{sample_id}.{library_id}.log",
    conda:
        "../../environments/coverm.yml"
    params:
        method=lambda w: w.method,
    shell:
        """
        ( coverm contig \
            --bam-files {input} \
            --methods {params.method} \
            --proper-pairs-only \
            --output-file /dev/stdout \
        | sed '1 s/Contig/sequence_id/' \
        | cut -f 1 -d " " \
        > {output} \
        ) 2> {log} 1>&2
        """


rule quantify__coverm__contig__join:
    """Aggregate coverm contig results"""
    input:
        lambda w: [
            COVERM
            / "contig"
            / f"{w.mag_catalogue}.{w.method}.{sample_id}.{library_id}.tsv"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        COVERM / "contig.{mag_catalogue}.{method}.tsv.gz",
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
