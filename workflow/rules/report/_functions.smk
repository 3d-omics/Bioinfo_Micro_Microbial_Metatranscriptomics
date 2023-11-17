def get_star_for_library_report(wildcards):
    """Get all star reports for a single library"""
    sample = wildcards.sample
    library = wildcards.library
    files = [STAR / f"{sample}.{library}.Log.final.out"]
    return files


def get_kraken2_for_library_report(wildcards):
    """Get all kraken2 reports for a single library"""
    sample = wildcards.sample
    library = wildcards.library
    files = [
        KRAKEN2 / kraken2_db / f"{sample}.{library}.report"
        for kraken2_db in KRAKEN2_DBS
    ]
    return files

def get_samtools_for_library_report(wildcards):
    sample = wildcards.sample
    library = wildcards.library
    files = [
        BOWTIE2 / f"{mag_catalogue}.{sample}.{library}.{report}"
        for mag_catalogue in MAG_CATALOGUES
        for report in BAM_REPORTS
    ]
    return files
