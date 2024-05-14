def get_mags_fasta(wildcards):
    """Get the path to the mag catalogue"""
    return features["mag_catalogues"][wildcards.mag_catalogue]["fasta"]


def get_mags_annotation(wildcards):
    """Get the path to the mag annotation in BED6 format"""
    return features["mag_catalogues"][wildcards.mag_catalogue]["bed6"]


def get_host_genome(wildcards):
    """Get the path to the host reference genome"""
    return features["hosts"][wildcards.host_name]["genome"]


def get_host_annotation(wildcards):
    """Get the path to the host annotation in GTF format"""
    return features["hosts"][wildcards.host_name]["gtf"]
