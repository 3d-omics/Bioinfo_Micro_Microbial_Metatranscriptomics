def compose_rg_id(wildcards):
    """Compose read group ID for bowtie2"""
    return f"{wildcards.sample}_{wildcards.library}"


def compose_rg_extra(wildcards):
    """Compose read group extra information for bowtie2"""
    return f"LB:truseq_{wildcards.library}\tPL:Illumina\tSM:{wildcards.sample}"


def get_method(wildcards):
    """Get the method wildcard"""
    return wildcards.method


def get_tsvs_for_assembly_coverm_genome(wildcards):
    """Get all the concrete coverm genome tsv files for a concrete assembly"""
    method = wildcards.method
    mag_catalogue = wildcards.mag_catalogue
    tsv_files = [
       COVERM / "genome" / mag_catalogue / method / f"{sample_id}.{library_id}.tsv"
        for sample_id, library_id in SAMPLE_LIB
    ]
    return tsv_files


def get_tsvs_for_assembly_coverm_contig(wildcards):
    """Get all the concrete coverm contig tsv files for a concrete assembly"""
    method = wildcards.method
    mag_catalogue = wildcards.mag_catalogue
    tsv_files = [
        COVERM / "contig" / mag_catalogue / method / f"{sample_id}.{library_id}.tsv"
        for sample_id, library_id in SAMPLE_LIB
    ]
    return tsv_files


def compose_input_dir_for_coverm_contig_aggregate(wildcards):
    mag_catalogue = wildcards.mag_catalogue
    method = wildcards.method
    return COVERM / "contig" / mag_catalogue / method


def compose_input_dir_for_coverm_genome_aggregate(wildcards):
    mag_catalogue = wildcards.mag_catalogue
    method = wildcards.method
    return COVERM / "genome" / mag_catalogue / method
