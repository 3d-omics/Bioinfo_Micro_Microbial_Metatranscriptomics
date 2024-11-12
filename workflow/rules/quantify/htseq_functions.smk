def get_tsvs_for_htseq(wildcards):
    mag_catalogue = wildcards.mag_catalogue
    tsv_files = [
        HTSEQ / mag_catalogue / f"{sample_id}.{library_id}.tsv.gz"
        for sample_id, library_id in SAMPLE_LIBRARY
    ]
    return tsv_files


def get_tsvs_for_subread(wildcards):
    mag_catalogue = wildcards.mag_catalogue
    tsv_files = [
        SUBREAD / mag_catalogue / f"{sample_id}.{library_id}.tsv.gz"
        for sample_id, library_id in SAMPLE_LIBRARY
    ]
    return tsv_files
