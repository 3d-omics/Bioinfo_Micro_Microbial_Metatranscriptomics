rule quantify__samtools__stats_cram__:
    """Compute stats for a cram"""
    input:
        cram=BOWTIE2 / "{mag_catalogue}.{sample}.{library}.cram",
        crai=BOWTIE2 / "{mag_catalogue}.{sample}.{library}.cram.crai",
        reference=REFERENCE / "mags" / "{mag_catalogue}.fa.gz",
        fai=REFERENCE / "mags" / "{mag_catalogue}.fa.gz.fai",
    output:
        tsv=BOWTIE2 / "{mag_catalogue}.{sample}.{library}.stats.tsv",
    log:
        BOWTIE2 / "{mag_catalogue}.{sample}.{library}.stats.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools stats --reference {input.reference} {input.cram} > {output.tsv} 2> {log}"
