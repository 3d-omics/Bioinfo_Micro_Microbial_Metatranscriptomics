def get_forward_for_bowtie2(wildcards):
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    if len(HOST_NAMES) == 0:
        return RIBODETECTOR / f"{sample_id}.{library_id}_1.fq.gz"
    return STAR / LAST_HOST / f"{sample_id}.{library_id}.Unmapped.out.mate1.gz"


def get_reverse_for_bowtie2(wildcards):
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    if len(HOST_NAMES) == 0:
        return RIBODETECTOR / f"{sample_id}.{library_id}_2.fq.gz"
    return STAR / LAST_HOST / f"{sample_id}.{library_id}.Unmapped.out.mate2.gz"


def compose_rg_id(wildcards):
    """Compose read group ID for bowtie2"""
    return f"{wildcards.sample_id}_{wildcards.library_id}"


def compose_rg_extra(wildcards):
    """Compose read group extra information for bowtie2"""
    return f"LB:truseq_{wildcards.library_id}\tPL:Illumina\tSM:{wildcards.sample_id}"


def get_method(wildcards):
    """Get the method wildcard"""
    return wildcards.method


def get_tsvs_for_assembly_coverm_genome(wildcards):
    """Get all the concrete coverm genome tsv files for a concrete assembly"""
    method = wildcards.method
    mag_catalogue = wildcards.mag_catalogue
    tsv_files = [
       COVERM / "genome" / mag_catalogue / method / f"{sample_id}.{library_id}.tsv"
        for sample_id, library_id in SAMPLE_LIBRARY
    ]
    return tsv_files


def get_tsvs_for_assembly_coverm_contig(wildcards):
    """Get all the concrete coverm contig tsv files for a concrete assembly"""
    method = wildcards.method
    mag_catalogue = wildcards.mag_catalogue
    tsv_files = [
        COVERM / "contig" / mag_catalogue / method / f"{sample_id}.{library_id}.tsv"
        for sample_id, library_id in SAMPLE_LIBRARY
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
