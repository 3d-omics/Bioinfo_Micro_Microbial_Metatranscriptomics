include: "__functions__.smk"
include: "mags.smk"


rule reference:
    """Prepare all the host and MAG references"""
    input:
        rules.reference__mags.input,
