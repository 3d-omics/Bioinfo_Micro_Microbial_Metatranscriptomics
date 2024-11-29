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
    return [STAR / f"{prev_genome}.{sample_id}.{library_id}_u{end}.fq.gz"]


def get_input_forward_for_host_mapping(wildcards):
    """Compose the forward input file"""
    return get_input_for_host_mapping(wildcards, "forward")


def get_input_reverse_for_host_mapping(wildcards):
    """Compose the forward input file"""
    return get_input_for_host_mapping(wildcards, "reverse")
