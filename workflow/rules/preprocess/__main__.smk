

include: "reads.smk"
include: "hosts.smk"
include: "fastp.smk"
include: "ribodetector.smk"
include: "kraken2.smk"
include: "star.smk"
include: "multiqc.smk"


rule preprocess__all:
    """Run the preprocessing steps:
    - fastp: trimming and adapter removal
    - ribodetector: removal of rRNAs
    - kraken2: screening of sequences
    - star: remove host RNA (if host is present)
    """
    input:
        rules.preprocess__fastp__all.input,
        rules.preprocess__ribodetector__all.input,
        rules.preprocess__kraken2__all.input,
        rules.preprocess__star.input,
        rules.preprocess__multiqc.output,
