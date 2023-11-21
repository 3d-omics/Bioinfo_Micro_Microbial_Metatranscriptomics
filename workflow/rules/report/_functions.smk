def get_star_for_library_report(wildcards):
    """Get all star reports for a single library"""
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    files = [
        STAR / host_name / f"{sample_id}.{library_id}.Log.final.out"
        for host_name in HOST_NAMES
    ]
    return files


def get_kraken2_for_library_report(wildcards):
    """Get all kraken2 reports for a single library"""
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    files = [
        KRAKEN2 / kraken2_db / f"{sample_id}.{library_id}.report"
        for kraken2_db in KRAKEN2_DBS
    ]
    return files


def get_samtools_for_library_report(wildcards):
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    files = [
        BOWTIE2 / f"{mag_catalogue}.{sample_id}.{library_id}.{report}"
        for mag_catalogue in MAG_CATALOGUES
        for report in BAM_REPORTS
    ]
    return files
