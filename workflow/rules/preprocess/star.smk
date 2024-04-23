

rule preprocess__star__align__:
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
        cram=STAR / "{host_name}" / "{sample_id}.{library_id}.cram",
        crai=STAR / "{host_name}" / "{sample_id}.{library_id}.cram.crai",
        u1=STAR / "{host_name}" / "{sample_id}.{library_id}.Unmapped.out.mate1.gz",
        u2=STAR / "{host_name}" / "{sample_id}.{library_id}.Unmapped.out.mate2.gz",
        report=STAR / "{host_name}" / "{sample_id}.{library_id}.Log.final.out",
        counts=STAR / "{host_name}" / "{sample_id}.{library_id}.ReadsPerGene.out.tab",
    log:
        STAR / "{host_name}" / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    params:
        out_prefix=get_star_out_prefix,
        u1=get_star_output_r1,
        u2=get_star_output_r2,
        bam=get_star_output_bam,
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

        pigz \
            --processes {threads} \
            --verbose \
            --best \
            {params.u1} \
            {params.u2} \
        2>> {log} 1>&2

        samtools view \
            --output-fmt CRAM \
            --reference {input.reference} \
            --threads {threads} \
            --write-index \
            --output {output.cram} \
            {params.bam} \
        2>> {log} 1>&2

        rm --verbose --force {params.bam} 2>> {log} 1>&2
        """


rule preprocess__star__align:
    """Get all the STAR counts for all hosts"""
    input:
        [
            STAR / host_name / f"{sample_id}.{library_id}.ReadsPerGene.out.tab"
            for sample_id, library_id in SAMPLE_LIBRARY
            for host_name in HOST_NAMES
        ],


rule preprocess__star__report:
    """Collect star reports"""
    input:
        [
            STAR / host_name / f"{sample_id}.{library_id}.Log.final.out"
            for sample_id, library_id in SAMPLE_LIBRARY
            for host_name in HOST_NAMES
        ],


rule preprocess__star:
    """Run all the elements in the star subworkflow"""
    input:
        rules.preprocess__star__align.input,
        rules.preprocess__star__report.input,
