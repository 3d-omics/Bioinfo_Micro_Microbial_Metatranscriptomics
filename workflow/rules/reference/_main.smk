include: "_functions.smk"
include: "hosts.smk"
include: "mags.smk"


rule reference:
    input:
        rules.reference__mags.input,
        rules.reference__hosts.input,
