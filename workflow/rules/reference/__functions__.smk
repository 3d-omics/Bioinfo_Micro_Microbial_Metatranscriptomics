def get_mags_fasta(wildcards):
    """Get the path to the mag catalogue"""
    return features["mag_catalogues"][wildcards.mag_catalogue]["fasta"]


def get_mags_bed6(wildcards):
    """Get the path to the mag annotation in BED6 format"""
    return features["mag_catalogues"][wildcards.mag_catalogue]["bed6"]


def get_mags_gtf(wildcards):
    """Get the path to the mag annotation in GTF format"""
    return features["mag_catalogues"][wildcards.mag_catalogue]["gtf"]
