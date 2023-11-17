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


def get_star_out_prefix(wildcards):
    """Get the star output folder from the library wildcards"""
    return STAR / f"{wildcards.sample}.{wildcards.library}."


def get_star_output_r1(wildcards):
    """Get the forward read output from the library wildcards"""
    return STAR / f"{wildcards.sample}.{wildcards.library}.Unmapped.out.mate1"


def get_star_output_r2(wildcards):
    """Get the reverse read output from the library wildcards"""
    return STAR / f"{wildcards.sample}.{wildcards.library}.Unmapped.out.mate2"


def get_kraken2_database(wildcards):
    """Get the path to the kraken2 database to be used"""
    return features["kraken2_databases"][wildcards.kraken2_db]
