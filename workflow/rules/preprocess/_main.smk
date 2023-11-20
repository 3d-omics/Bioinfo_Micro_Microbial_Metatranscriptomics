include: "_functions.smk"
include: "fastp.smk"
include: "ribodetector.smk"
include: "kraken2.smk"
include: "star.smk"


rule preprocessing:
    input:
        rules.fastp.input,
        rules.ribodetector.input,
        rules.kraken2.input,
        rules.star.input,
