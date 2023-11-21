include: "_functions.smk"
include: "hosts.smk"
include: "mags.smk"


rule reference:
    """Prepare all the host and MAG references"""
    input:
        rules.reference__mags.input,
        rules.reference__hosts.input,
