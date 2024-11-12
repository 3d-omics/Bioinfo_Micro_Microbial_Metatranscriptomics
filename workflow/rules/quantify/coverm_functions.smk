def get_method(wildcards):
    """Get the method wildcard"""
    return wildcards.method


def get_min_covered_fraction(wildcards):
    """Get the min covered fraction"""
    return params["quantify"]["coverm"]["genome"]["min_covered_fraction"]


def get_separator(wildcards):
    """Get the mag-contig separator"""
    return params["quantify"]["coverm"]["separator"]


def get_tsvs_for_assembly_coverm(wildcards, genome_or_contig):
    """Get all the concrete coverm tsv files for a concrete assembly"""
    assert genome_or_contig in ["genome", "contig"]
    method = wildcards.method
    mag_catalogue = wildcards.mag_catalogue
    tsv_files = [
        COVERM
        / genome_or_contig
        / mag_catalogue
        / method
        / f"{sample_id}.{library_id}.tsv.gz"
        for sample_id, library_id in SAMPLE_LIBRARY
    ]
    return tsv_files


def get_tsvs_for_assembly_coverm_genome(wildcards):
    """Get all the concrete coverm genome tsv files for a concrete assembly"""
    return get_tsvs_for_assembly_coverm(wildcards, "genome")


def get_tsvs_for_assembly_coverm_contig(wildcards):
    """Get all the concrete coverm contig tsv files for a concrete assembly"""
    return get_tsvs_for_assembly_coverm(wildcards, "contig")


def compose_input_dir_for_coverm_contig_aggregate(wildcards):
    return COVERM / "contig" / wildcards.mag_catalogue / wildcards.method


def compose_input_dir_for_coverm_genome_aggregate(wildcards):
    return COVERM / "genome" / wildcards.mag_catalogue / wildcards.method
