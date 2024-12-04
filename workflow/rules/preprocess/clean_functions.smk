def get_final_fastq(wildcards, forward_or_reverse):
    """Get the forward or reverse for bowtie2"""
    assert forward_or_reverse in ["forward", "reverse"]
    end = 1 if forward_or_reverse == "forward" else 2
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    if len(HOST_NAMES) == 0:
        return PRE_RIBODETECTOR / f"{sample_id}.{library_id}_{end}.fq.gz"
    return PRE_STAR / LAST_HOST / f"{sample_id}.{library_id}_u{end}.fq.gz"


def get_final_forward(wildcards):
    """Get the forward for bowtie2"""
    return get_final_fastq(wildcards, "forward")


def get_final_reverse(wildcards):
    """Get the reverse for bowtie2"""
    return get_final_fastq(wildcards, "reverse")
