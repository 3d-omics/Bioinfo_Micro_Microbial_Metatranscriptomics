include: "star_functions.smk"


rule preprocess__star__index:
    """Index the genome for STAR"""
    input:
        genome=HOSTS / "{host_name}.fa",
        annotation=HOSTS / "{host_name}.gtf",
    output:
        multiext(
            str(STAR_INDEX / "{host_name}") + "/",
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
        "../../environments/star.yml"
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


rule preprocess__star__index__all:
    """Build all the STAR indexes"""
    input:
        [STAR_INDEX / host_name / "Genome" for host_name in HOST_NAMES],


rule preprocess__star__align:
    """Align one library to the host genome with STAR to discard host RNA"""
    input:
        forward_=get_input_forward_for_host_mapping,
        reverse_=get_input_reverse_for_host_mapping,
        index=multiext(
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
        reference=HOSTS / "{host_name}.fa",
        fai=HOSTS / "{host_name}.fa.fai",
    output:
        bam=STAR / "{host_name}.{sample_id}.{library_id}.Aligned.sortedByCoord.out.bam",
        report=STAR / "{host_name}.{sample_id}.{library_id}.Log.final.out",
        counts=STAR / "{host_name}.{sample_id}.{library_id}.ReadsPerGene.out.tab",
    log:
        STAR / "{host_name}.{sample_id}.{library_id}.log",
    conda:
        "../../environments/star.yml"
    params:
        out_prefix=get_star_out_prefix,
        index=lambda w: STAR_INDEX / w.host_name,
    retries: 5
    shell:
        """
        ulimit -n 90000 2> {log} 1>&2

        STAR \
            --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {params.index} \
            --readFilesIn \
                {input.forward_} \
                {input.reverse_} \
            --outFileNamePrefix {params.out_prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within KeepPairs \
            --outReadsUnmapped Fastx \
            --readFilesCommand "gzip -cd" \
            --quantMode GeneCounts \
        2>> {log} 1>&2
        """


rule preprocess__star__align__all:
    """Get all the STAR counts for all hosts"""
    input:
        [
            STAR
            / f"{host_name}.{sample_id}.{library_id}.Aligned.sortedByCoord.out.bam"
            for sample_id, library_id in SAMPLE_LIBRARY
            for host_name in HOST_NAMES
        ],


rule preprocess__star__fastq:
    """Convert BAM to FASTQ using samtools and using the correct reference

    NOTE: bowtie2 does not like CRAM files, and although can use a BAM file as an input,
    bowtie2 fails to receive a piped SAM input. Therefore, we need to convert the CRAM file to a physical FASTQ file.
    """
    input:
        bam=STAR / "{host_name}.{sample_id}.{library_id}.Aligned.sortedByCoord.out.bam",
        bai=STAR
        / "{host_name}.{sample_id}.{library_id}.Aligned.sortedByCoord.out.bam.bai",
    output:
        forward_=temp(STAR / "{host_name}.{sample_id}.{library_id}_u1.fq.gz"),
        reverse_=temp(STAR / "{host_name}.{sample_id}.{library_id}_u2.fq.gz"),
    log:
        STAR / "{host_name}.{sample_id}.{library_id}.unaligned.log",
    conda:
        "../../environments/star.yml"
    shell:
        """
        rm \
            --recursive \
            --force \
            {output.forward_}.collate

        ( samtools view \
            -f 12 \
            -u \
            --threads {threads} \
            {input} \
            "*" \
        | samtools collate \
            -O \
            -u \
            -f \
            -r 1e6 \
            -T {output.forward_}.collate \
            --threads {threads} \
            - \
        | samtools fastq \
            -1 {output.forward_} \
            -2 {output.reverse_} \
            -0 /dev/null \
            -s /dev/null \
            --threads {threads} \
            -c 1 \
            /dev/stdin \
        ) 2> {log} 1>&2
        """


rule preprocess__star__all:
    """Run all the elements in the star subworkflow"""
    input:
        rules.preprocess__star__index__all.input,
        rules.preprocess__star__align__all.input,
