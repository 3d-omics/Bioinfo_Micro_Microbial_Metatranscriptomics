# bowtie2 ----
def get_file_for_bowtie2(wildcards, forward_or_reverse):
    """Get the forward or reverse for bowtie2"""
    assert forward_or_reverse in ["forward", "reverse"]
    end = 1 if forward_or_reverse == "forward" else 2
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    if len(HOST_NAMES) == 0:
        return RIBODETECTOR / f"{sample_id}.{library_id}_{end}.fq.gz"
    return STAR / LAST_HOST / f"{sample_id}.{library_id}.Unmapped.out.mate{end}.gz"


def get_forward_for_bowtie2(wildcards):
    """Get the forward for bowtie2"""
    return get_file_for_bowtie2(wildcards, "forward")


def get_reverse_for_bowtie2(wildcards):
    """Get the reverse for bowtie2"""
    return get_file_for_bowtie2(wildcards, "reverse")


def compose_rg_id(wildcards):
    """Compose read group ID for bowtie2"""
    return f"{wildcards.sample_id}_{wildcards.library_id}"


def compose_rg_extra(wildcards):
    """Compose read group extra information for bowtie2"""
    return f"LB:truseq_{wildcards.library_id}\tPL:Illumina\tSM:{wildcards.sample_id}"


# coverm ----
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
        / f"{sample_id}.{library_id}.tsv"
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


# Bedtools
def get_tsvs_for_bedtools(wildcards):
    mag_catalogue = wildcards.mag_catalogue
    tsv_files = [
        BEDTOOLS / mag_catalogue / f"{sample_id}.{library_id}.tsv"
        for sample_id, library_id in SAMPLE_LIBRARY
    ]
    return tsv_files


# htseq
def get_tsvs_for_htseq(wildcards):
    mag_catalogue = wildcards.mag_catalogue
    tsv_files = [
        HTSEQ / mag_catalogue / f"{sample_id}.{library_id}.tsv"
        for sample_id, library_id in SAMPLE_LIBRARY
    ]
    return tsv_files


def get_tsvs_for_subread(wildcards):
    mag_catalogue = wildcards.mag_catalogue
    tsv_files = [
        SUBREAD / mag_catalogue / f"{sample_id}.{library_id}.tsv"
        for sample_id, library_id in SAMPLE_LIBRARY
    ]
    return tsv_files
