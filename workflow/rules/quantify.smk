include: "quantify/mags.smk"
include: "quantify/bowtie2.smk"
include: "quantify/subread.smk"
include: "quantify/multiqc.smk"


rule quantify__all:
    """
    Run the quantify steps:
    - bowtie2: Map to MAG catalogue
    - subread: Quantify reads per CDS, tRNA and rRNA in the GTF file
    """
    input:
        rules.quantify__mags__all.input,
        rules.quantify__bowtie2__all.input,
        rules.quantify__subread__all.input,
        rules.quantify__multiqc__all.input,
