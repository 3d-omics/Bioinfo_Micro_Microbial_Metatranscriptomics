def get_host_genome(wildcards):
    """Get the path to the host reference genome"""
    return features["hosts"][wildcards.host_name]["genome"]


def get_host_annotation(wildcards):
    """Get the path to the host annotation in GTF format"""
    return features["hosts"][wildcards.host_name]["gtf"]
