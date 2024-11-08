include: "__functions__.smk"
include: "mags.smk"
include: "bowtie2.smk"
include: "coverm.smk"
include: "samtools.smk"
# include: "bedtools.smk"
include: "htseq.smk"
include: "subread.smk"
include: "multiqc.smk"


rule quantify:
    """
    Run the quantify steps:
    - bowtie2: Map to MAG catalogue
    - coverm: Get the count tables
    """
    input:
        rules.quantify__mags__all.input,
        rules.quantify__bowtie2__all.input,
        rules.quantify__coverm.input,
        # rules.quantify__bedtools.input,
        rules.quantify__htseq.input,
        rules.quantify__subread.input,
        rules.quantify__multiqc__all.input,
