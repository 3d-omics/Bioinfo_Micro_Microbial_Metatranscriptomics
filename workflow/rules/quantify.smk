include: "quantify/mags.smk"
include: "quantify/bowtie2.smk"
include: "quantify/coverm.smk"
include: "quantify/samtools.smk"
include: "quantify/htseq.smk"
include: "quantify/subread.smk"
include: "quantify/multiqc.smk"


rule quantify__all:
    """
    Run the quantify steps:
    - bowtie2: Map to MAG catalogue
    - coverm: Get the count tables
    """
    input:
        rules.quantify__mags__all.input,
        rules.quantify__bowtie2__all.input,
        rules.quantify__coverm__all.input,
        rules.quantify__htseq__all.input,
        rules.quantify__subread__all.input,
        rules.quantify__multiqc__all.input,
