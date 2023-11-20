rule kraken2_assign_one:
    """Run kraken2 over one library

    The database must be provided by the user in the config file.
    """
    input:
        forward_=FASTP / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=FASTP / "{sample_id}.{library_id}_2.fq.gz",
        database=get_kraken2_database,
    output:
        out_gz=KRAKEN2 / "{kraken2_db}" / "{sample_id}.{library_id}.out.gz",
        report=KRAKEN2 / "{kraken2_db}" / "{sample_id}.{library_id}.report",
    log:
        log=KRAKEN2 / "{kraken2_db}" / "{sample_id}.{library_id}.log",
    conda:
        "_env.yml"
    threads: 24
    resources:
        mem_mb=eval(params["kraken2"]["mem_mb"]),
        runtime=24 * 60,
    shell:
        """
        kraken2 \
            --db {input.database} \
            --threads {threads} \
            --paired \
            --gzip-compressed \
            --output >(pigz -11 > {output.out_gz}) \
            --report {output.report} \
            {input.forward_} \
            {input.reverse_} \
        > {log} 2>&1
        """


rule kraken2_assign_all:
    """Run kraken2 over all libraries"""
    input:
        [
            KRAKEN2 / f"{kraken2_db}" / f"{sample_id}.{library_id}.report"
            for sample_id, library_id in SAMPLE_LIBRARY
            for kraken2_db in KRAKEN2_DBS
        ],


rule kraken2_report_one:
    """Generate a report for one library

    Equivalent to just runing kraken2.
    """
    input:
        rules.kraken2_assign_one.output.report,


rule kraken2_report_all:
    """Generate all the reports"""
    input:
        rules.kraken2_assign_all.input,


rule kraken2:
    """Run the kraken2 subworkflow"""
    input:
        rules.kraken2_assign_all.input,
        rules.kraken2_report_all.input,
