def get_forward_adapter(wildcards):
    """Get forward adapter for a sample and library."""
    forward_adapter =  samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["forward_adapter"].tolist()[0]
    if pd.isna(forward_adapter):
        return "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    return forward_adapter



def get_reverse_adapter(wildcards):
    """Get reverse adapter for a sample and library."""
    reverse_adapter = samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["reverse_adapter"].tolist()[0]
    if pd.isna(reverse_adapter):
        return "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    return reverse_adapter
