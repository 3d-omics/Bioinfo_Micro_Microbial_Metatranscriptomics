include: "reads_functions.smk"


rule preprocess__reads:
    """Make a link to the original file, with a prettier name than default"""
    input:
        forward_=get_forward,
        reverse_=get_reverse,
    output:
        forward_=PRE_READS / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_READS / "{sample_id}.{library_id}_2.fq.gz",
    log:
        PRE_READS / "{sample_id}.{library_id}.log",
    conda:
        "base"
    group:
        "preprocess__{sample_id}.{library_id}"
    shell:
        """
        ln --symbolic $(readlink --canonicalize {input.forward_}) {output.forward_}
        ln --symbolic $(readlink --canonicalize {input.reverse_}) {output.reverse_}
        """


rule preprocess__reads__link__all:
    """Link all reads in the samples.tsv"""
    input:
        [
            PRE_READS / f"{sample_id}.{library_id}_{end}.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in ["1", "2"]
        ],


rule preprocess__reads__fastqc__all:
    """Run fastqc on all raw reads"""
    input:
        [
            PRE_READS / f"{sample_id}.{library_id}_{end}_fastqc.{extension}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in ["1", "2"]
            for extension in ["html", "zip"]
        ],


rule preprocess__reads__all:
    input:
        rules.preprocess__reads__link__all.input,
        rules.preprocess__reads__fastqc__all.input,
