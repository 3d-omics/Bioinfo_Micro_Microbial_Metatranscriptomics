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
