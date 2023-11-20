rule star_index_one:
    """Index the genome for STAR"""
    input:
        genome=HOSTS / "{host_name}.fa",
        annotation=HOSTS / "{host_name}.gtf",
    output:
        folder=directory(STAR_INDEX / "{host_name}"),
    params:
        sjdbOverhang=params["star"]["index"]["sjdbOverhang"],
    conda:
        "_env.yml"
    log:
        STAR_INDEX / "{host_name}.log"
    threads: 24
    resources:
        mem_mb=32 * 1024,
        runtime=24 * 60,
    shell:
        """
        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output.folder} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.annotation} \
            --sjdbOverhang {params.sjdbOverhang} \
        2> {log} 1>&2
        """


rule star_index:
    input:
        [
            STAR_INDEX / f"{host_name}"
            for host_name in HOST_NAMES
        ]


rule star_align_one:
    """Align one library to the host genome with STAR to discard host RNA"""
    input:
        forward_=get_input_forward_for_host_mapping,
        reverse_=get_input_reverse_for_host_mapping,
        index=STAR_INDEX / "{host_name}",
    output:
        bam=temp(STAR / "{host_name}" / "{sample_id}.{library_id}.Aligned.sortedByCoord.out.bam"),
        u1=temp(STAR / "{host_name}" / "{sample_id}.{library_id}.Unmapped.out.mate1"),
        u2=temp(STAR / "{host_name}" / "{sample_id}.{library_id}.Unmapped.out.mate2"),
        report=STAR / "{host_name}" / "{sample_id}.{library_id}.Log.final.out",
        counts=STAR / "{host_name}" / "{sample_id}.{library_id}.ReadsPerGene.out.tab",
    log:
        STAR / "{host_name}" / "{sample_id}.{library_id}.log",
    params:
        out_prefix=get_star_out_prefix,
    conda:
        "_env.yml"
    threads: 24
    resources:
        mem_mb = 32 * 1024,
        runtime=24 * 60,
    shell:
        """
        ulimit -n 90000 2> {log} 1>&2

        STAR \
            --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {input.index} \
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


rule star_align_all:
    input:
        [
            STAR / host_name / f"{sample_id}.{library_id}.ReadsPerGene.out.tab"
            for sample_id, library_id in SAMPLE_LIBRARY
            for host_name in HOST_NAMES
        ],


rule star_cram_one:
    """Convert to cram one library

    Note: we use samtools sort when it is already sorted because there is no
    other way to use minimizers on the unmapped fraction.
    """
    input:
        bam=STAR / "{host_name}" / "{sample_id}.{library_id}.Aligned.sortedByCoord.out.bam",
        reference=HOSTS / "{host_name}.fa",
        fai = HOSTS / "{host_name}.fa.fai",
    output:
        cram=STAR / "{host_name}" / "{sample_id}.{library_id}.cram",
        crai=STAR / "{host_name}" / "{sample_id}.{library_id}.cram.crai",
    log:
        STAR / "{host_name}" / "{sample_id}.{library_id}.Aligned.sortedByCoord.out.cram.log",
    conda:
        "_env.yml"
    threads: 24
    resources:
        mem_mb=24 * 1024,
        runtime=24 * 60,
    shell:
        """
        samtools sort \
            -l 9 \
            -m 1G \
            -M \
            -o {output.cram} \
            --output-fmt CRAM \
            --reference {input.reference} \
            -@ {threads} \
            --write-index \
            {input.bam} \
        2> {log} 1>&2
        """


rule star_cram_all:
    """Convert to cram all libraries"""
    input:
        [
            STAR / host_name / f"{sample_id}.{library_id}.cram"
            for sample_id, library_id in SAMPLE_LIBRARY
            for host_name in HOST_NAMES
        ],


rule star_report_all:
    """Collect star reports"""
    input:
        [
            STAR / host_name / f"{sample_id}.{library_id}.Log.final.out"
            for sample_id, library_id in SAMPLE_LIBRARY
            for host_name in HOST_NAMES
        ],


rule star:
    """Run all the elements in the star subworkflow"""
    input:
        rules.star_align_all.input,
        rules.star_cram_all.input,
        rules.star_report_all.input,
