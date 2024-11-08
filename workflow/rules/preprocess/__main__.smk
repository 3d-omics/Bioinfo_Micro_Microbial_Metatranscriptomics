include: "reads_functions.smk"
include: "reads.smk"
include: "__functions__.smk"
include: "fastp.smk"
include: "ribodetector.smk"
include: "kraken2.smk"
include: "index.smk"
include: "star.smk"


rule preprocess:
    """Run the preprocessing steps:
    - fastp: trimming and adapter removal
    - ribodetector: removal of rRNAs
    - kraken2: screening of sequences
    - star: remove host RNA (if host is present)
    """
    input:
        rules.preprocess__fastp.input,
        rules.preprocess__ribodetector.input,
        rules.preprocess__kraken2.input,
        rules.preprocess__star.input,
