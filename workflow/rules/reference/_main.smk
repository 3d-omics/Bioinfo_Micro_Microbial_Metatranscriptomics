include: "_functions.smk"
include: "hosts.smk"
include: "mags.smk"


rule reference:
    input:
        rules.reference_mags.input,
        rules.reference_hosts.input,
