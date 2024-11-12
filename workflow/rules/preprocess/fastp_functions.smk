# adapters ----
def get_adapter(wildcards, forward_or_reverse):
    """Get the correct adapter"""
    assert forward_or_reverse in ["forward", "reverse"]
    end = "forward_adapter" if forward_or_reverse == "forward" else "reverse_adapter"
    return samples[
        (samples["sample_id"] == wildcards.sample_id)
        & (samples["library_id"] == wildcards.library_id)
    ][end].tolist()[0]


def get_forward_adapter(wildcards):
    """Get forward adapter for a sample and library."""
    return get_adapter(wildcards, "forward")


def get_reverse_adapter(wildcards):
    """Get reverse adapter for a sample and library."""
    return get_adapter(wildcards, "reverse")


def compose_adapters(wildcards):
    """Compose the forward and reverse adapter line for fastp"""
    forward = get_forward_adapter(wildcards)
    reverse = get_reverse_adapter(wildcards)
    return f"--adapter_sequence {forward} --adapter_sequence_r2 {reverse}"
