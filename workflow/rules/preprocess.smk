

include: "preprocess/reads.smk"
include: "preprocess/hosts.smk"
include: "preprocess/fastp.smk"
include: "preprocess/ribodetector.smk"
include: "preprocess/kraken2.smk"
include: "preprocess/star.smk"
include: "preprocess/multiqc.smk"


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
        rules.preprocess__star__all.input,
        rules.preprocess__multiqc__all.output,
