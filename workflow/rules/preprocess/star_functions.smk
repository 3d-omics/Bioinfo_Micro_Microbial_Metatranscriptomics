def get_input_for_host_mapping(wildcards, forward_or_reverse):
    """Get the forward or reverse file for host mapping"""
    assert forward_or_reverse in ["forward", "reverse"]
    end = 1 if forward_or_reverse == "forward" else 2
    host_name = wildcards.host_name
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    if host_name == HOST_NAMES[0]:
        return RIBODETECTOR / f"{sample_id}.{library_id}_{end}.fq.gz"
    genome_index = HOST_NAMES.index(host_name)
    prev_genome = HOST_NAMES[genome_index - 1]
    return [STAR / prev_genome / f"{sample_id}.{library_id}.Unmapped.out.mate{end}.gz"]


def get_input_forward_for_host_mapping(wildcards):
    """Compose the forward input file"""
    return get_input_for_host_mapping(wildcards, "forward")


def get_input_reverse_for_host_mapping(wildcards):
    """Compose the forward input file"""
    return get_input_for_host_mapping(wildcards, "reverse")


def get_star_out_prefix(wildcards):
    """Get the star output folder from the library wildcards"""
    return STAR / wildcards.host_name / f"{wildcards.sample_id}.{wildcards.library_id}."


def get_star_output_r(wildcards, forward_or_reverse):
    assert forward_or_reverse in ["forward", "reverse"]
    end = 1 if forward_or_reverse == "forward" else 2
    host_name = wildcards.host_name
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    return STAR / host_name / f"{sample_id}.{library_id}.Unmapped.out.mate{end}"


def get_star_output_r1(wildcards):
    """Get the forward read output from the library wildcards"""
    return get_star_output_r(wildcards, "forward")


def get_star_output_r2(wildcards):
    """Get the reverse read output from the library wildcards"""
    return get_star_output_r(wildcards, "reverse")


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
    suffix = "Aligned.sortedByCoord.out.bam"
    return STAR / host_name / f"{sample_id}.{library_id}.{suffix}"
