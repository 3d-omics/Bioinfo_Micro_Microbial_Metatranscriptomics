include: "_functions.smk"
include: "bowtie2.smk"
include: "coverm.smk"
include: "samtools.smk"


rule quantification:
    input:
        rules.coverm.input,
