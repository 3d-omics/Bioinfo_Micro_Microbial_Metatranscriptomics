include: "_functions.smk"
include: "bowtie2.smk"
include: "coverm.smk"
include: "samtools.smk"


rule quantify:
    input:
        rules.quantify__coverm.input,
