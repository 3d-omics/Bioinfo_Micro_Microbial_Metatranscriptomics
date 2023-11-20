include: "_functions.smk"
include: "fastp.smk"
include: "ribodetector.smk"
include: "kraken2.smk"
include: "star.smk"


rule preprocess:
    input:
        rules.preprocess__fastp.input,
        rules.preprocess__ribodetector.input,
        rules.preprocess__kraken2.input,
        rules.preprocess__star.input,
