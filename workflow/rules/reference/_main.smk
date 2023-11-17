include: "_functions.smk"
include: "hosts.smk"
include: "mags.smk"


rule reference:
    input:
        rules.reference_mags.input,
        rules.reference_set_dna.output.fa,
        rules.reference_set_gtf.output.gtf,
