# fastp ----
def get_adapter(wildcards, forward_or_reverse):
    """Get forward or reverse adapter"""
    assert forward_or_reverse in ["forward_adapter", "reverse_adapter"]
    return samples[
        (samples["sample_id"] == wildcards.sample_id)
        & (samples["library_id"] == wildcards.library_id)
    ][forward_or_reverse].tolist()[0]


def get_forward_adapter(wildcards):
    """Get forward adapter for a sample and library."""
    return get_adapter(wildcards, "forward_adapter")


def get_reverse_adapter(wildcards):
    """Get reverse adapter for a sample and library."""
    return get_adapter(wildcards, "reverse_adapter")


#
def get_input_forward_for_host_mapping(wildcards):
    """Compose the forward input file"""
    host_name = wildcards.host_name
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    if host_name == HOST_NAMES[0]:
        return RIBODETECTOR / f"{sample_id}.{library_id}_1.fq.gz"
    genome_index = HOST_NAMES.index(host_name)
    prev_genome = HOST_NAMES[genome_index - 1]
    return [STAR / prev_genome / f"{sample_id}.{library_id}.Unmapped.out.mate1.gz"]


def get_input_reverse_for_host_mapping(wildcards):
    """Compose the forward input file"""
    host_name = wildcards.host_name
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    if host_name == HOST_NAMES[0]:
        return RIBODETECTOR / f"{sample_id}.{library_id}_1.fq.gz"
    genome_index = HOST_NAMES.index(host_name)
    prev_genome = HOST_NAMES[genome_index - 1]
    return [STAR / prev_genome / f"{sample_id}.{library_id}.Unmapped.out.mate2.gz"]


# star ----
def get_star_out_prefix(wildcards):
    """Get the star output folder from the library wildcards"""
    return STAR / wildcards.host_name / f"{wildcards.sample_id}.{wildcards.library_id}."


def get_star_output_r1(wildcards):
    """Get the forward read output from the library wildcards"""
    return (
        STAR
        / wildcards.host_name
        / f"{wildcards.sample_id}.{wildcards.library_id}.Unmapped.out.mate1"
    )


def get_star_output_r2(wildcards):
    """Get the reverse read output from the library wildcards"""
    return (
        STAR
        / wildcards.host_name
        / f"{wildcards.sample_id}.{wildcards.library_id}.Unmapped.out.mate2"
    )


def get_star_output_r1_gz(wildcards):
    """get_star_output_r1 with gz"""
    return get_star_output_r1(wildcards) + ".gz"


def get_star_output_r2_gz(wildcards):
    """get_star_output_r2 with gz"""
    return get_star_output_r2(wildcards) + ".gz"


def get_star_output_bam(wildcards):
    """Get the star generated bam"""
    host_name = wildcards.host_name
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    return (
        STAR
        / f"{host_name}"
        / f"{sample_id}.{library_id}.Aligned.sortedByCoord.out.bam"
    )


# kraken2 ----
def get_kraken2_database(wildcards):
    """Get the path to the kraken2 database to be used"""
    return features["databases"]["kraken2"][wildcards.kraken2_db]


def compose_out_folder_for_eval_kraken2_assign_all(wildcards):
    """Just compose the output folder"""
    return KRAKEN2 / wildcards.kraken2_db
