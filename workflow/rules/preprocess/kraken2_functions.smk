def get_kraken2_database(wildcards):
    """Get the path to the kraken2 database to be used"""
    return features["databases"]["kraken2"][wildcards.kraken2_db]


def compose_out_folder_for_eval_kraken2_assign_all(wildcards):
    """Just compose the output folder"""
    return KRAKEN2 / wildcards.kraken2_db
