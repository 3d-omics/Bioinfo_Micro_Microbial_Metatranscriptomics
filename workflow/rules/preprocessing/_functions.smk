def get_forward_adapter(wildcards):
    """Get forward adapter for a sample and library."""
    forward_adapter =  samples[
        (samples["sample_id"] == wildcards.sample_id)
        & (samples["library_id"] == wildcards.library_id)
    ]["forward_adapter"].tolist()[0]
    if pd.isna(forward_adapter):
        return "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    return forward_adapter


def get_reverse_adapter(wildcards):
    """Get reverse adapter for a sample and library."""
    reverse_adapter = samples[
        (samples["sample_id"] == wildcards.sample_id)
        & (samples["library_id"] == wildcards.library_id)
    ]["reverse_adapter"].tolist()[0]
    if pd.isna(reverse_adapter):
        return "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    return reverse_adapter

def get_input_forward_for_host_mapping(wildcards):
    """Compose the forward input file"""
    host_name = wildcards.host_name
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    if host_name == HOST_NAMES[0]:
        return RIBODETECTOR / f"{sample_id}.{library_id}_1.fq.gz"
    genome_index = HOST_NAMES.index(host_name)
    prev_genome = HOST_NAMES[genome_index - 1]
    return [
        STAR / prev_genome / f"{sample_id}.{library_id}.Unmapped.out.mate1.gz"
    ]


def get_input_reverse_for_host_mapping(wildcards):
    """Compose the forward input file"""
    host_name = wildcards.host_name
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    if host_name == HOST_NAMES[0]:
        return RIBODETECTOR / f"{sample_id}.{library_id}_1.fq.gz"
    genome_index = HOST_NAMES.index(host_name)
    prev_genome = HOST_NAMES[genome_index - 1]
    return [
        STAR / prev_genome / f"{sample_id}.{library_id}.Unmapped.out.mate2.gz"
    ]


def get_star_out_prefix(wildcards):
    """Get the star output folder from the library wildcards"""
    return STAR / wildcards.host_name / f"{wildcards.sample_id}.{wildcards.library_id}."


def get_star_output_r1(wildcards):
    """Get the forward read output from the library wildcards"""
    return STAR /  wildcards.host_name / f"{wildcards.sample_id}.{wildcards.library_id}.Unmapped.out.mate1"


def get_star_output_r2(wildcards):
    """Get the reverse read output from the library wildcards"""
    return STAR /  wildcards.host_name / f"{wildcards.sample_id}.{wildcards.library_id}.Unmapped.out.mate2"

def get_star_output_r1_gz(wildcards):
    """get_star_output_r1 with gz"""
    return get_star_output_r1(wildcards) + ".gz"


def get_star_output_r2_gz(wildcards):
    """get_star_output_r2 with gz"""
    return  get_star_output_r2(wildcards) + ".gz"


def get_kraken2_database(wildcards):
    """Get the path to the kraken2 database to be used"""
    return features["kraken2_databases"][wildcards.kraken2_db]
