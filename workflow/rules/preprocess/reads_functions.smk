def get_read_file(wildcards, forward_or_reverse):
    """Get the correct read file"""
    assert forward_or_reverse in ["forward", "reverse"]
    end = "forward_filename" if forward_or_reverse == "forward" else "reverse_filename"
    return samples[
        (samples["sample_id"] == wildcards.sample_id)
        & (samples["library_id"] == wildcards.library_id)
    ][end].tolist()[0]


def get_forward(wildcards):
    """Get forward reads for a sample and library."""
    return get_read_file(wildcards, "forward")


def get_reverse(wildcards):
    """Get reverse reads for a sample and library."""
    return get_read_file(wildcards, "reverse")
