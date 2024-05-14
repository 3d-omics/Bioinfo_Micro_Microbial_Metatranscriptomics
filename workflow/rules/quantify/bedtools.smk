rule quantify__bedtools__:
    input:
        cram=BOWTIE2 / "{mag_catalogue}.{sample_id}.{library_id}.cram",
        crai=BOWTIE2 / "{mag_catalogue}.{sample_id}.{library_id}.cram.crai",
        reference=MAGS / "{mag_catalogue}.fa.gz",
        fai=MAGS / "{mag_catalogue}.fa.gz.fai",
        annotation=MAGS / "{mag_catalogue}.bed6",
    output:
        BEDTOOLS / "{mag_catalogue}" / "{sample_id}.{library_id}.tsv",
    log:
        BEDTOOLS / "{mag_catalogue}" / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        echo -E "contig\tcount" > {output} 2> {log}

        ( samtools view \
            --uncompressed \
            --reference {input.reference} \
            {input.cram} \
        | bedtools intersect \
            -a {input.annotation} \
            -b - \
            -c \
            -bed \
        | cut -f 4,6 \
        | sort -k1,1 \
        | uniq --count \
        | awk -v OFS="\t" '{{print $2, $1}}' \
        >> {output} \
        ) 2>> {log}
        """


rule quantify__bedtools__aggregate__:
    input:
        tsvs=get_tsvs_for_bedtools,
    output:
        BEDTOOLS / "{mag_catalogue}.tsv",
    log:
        BEDTOOLS / "{mag_catalogue}.log",
    conda:
        "__environment__.yml"
    params:
        input_folder=lambda w: BEDTOOLS / f"{w.mag_catalogue}",
    shell:
        """
        Rscript --vanilla \
            workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_folder} \
            --output-file {output} \
        2> {log}
        """


rule quantify__bedtools:
    input:
        [BEDTOOLS / f"{mag_catalogue}.tsv" for mag_catalogue in MAG_CATALOGUES],
