rule preprocess__index__:
    """Index the genome for STAR"""
    input:
        genome=HOSTS / "{host_name}.fa",
        annotation=HOSTS / "{host_name}.gtf",
    output:
        multiext(
            str(STAR_INDEX) + "/{host_name}/",
            "chrLength.txt",
            "chrNameLength.txt",
            "chrName.txt",
            "chrStart.txt",
            "exonGeTrInfo.tab",
            "exonInfo.tab",
            "geneInfo.tab",
            "Genome",
            "genomeParameters.txt",
            "Log.out",
            "SA",
            "SAindex",
            "sjdbInfo.txt",
            "sjdbList.fromGTF.out.tab",
            "sjdbList.out.tab",
            "transcriptInfo.tab",
        ),
    conda:
        "__environment__.yml"
    log:
        STAR_INDEX / "{host_name}.log",
    params:
        sjdbOverhang=params["preprocess"]["star"]["index"]["sjdbOverhang"],
        prefix=lambda w: str(STAR_INDEX / w.host_name),
    retries: 5
    cache: True
    shell:
        """
        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {params.prefix} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.annotation} \
            --sjdbOverhang {params.sjdbOverhang} \
        2> {log} 1>&2
        """


rule preprocess__index:
    """Build all the STAR indexes"""
    input:
        [STAR_INDEX / host_name / "Genome" for host_name in HOST_NAMES],
