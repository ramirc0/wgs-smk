# CRAM conversion and its validation.


rule convert_to_cram:
    input:
        **REFERENCE,
        bam=f"{OUTDIR}/{{sample}}/sort_bam/{{sample}}.bam",
    output:
        cram=f"{OUTDIR}/{{sample}}/convert_to_cram/{{sample}}.cram",
        crai=f"{OUTDIR}/{{sample}}/convert_to_cram/{{sample}}.cram.crai",
        cram_md5=f"{OUTDIR}/{{sample}}/convert_to_cram/{{sample}}.cram.md5",
    log:
        f"{LOGDIR}/convert_to_cram/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/convert_to_cram/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        samtools view \
            -@ {threads:q} \
            --output-fmt-option version=3.0 \
            -C \
            -T {input.fasta:q} \
            --write-index \
            -o {output.cram:q} \
            {input.bam:q}
        md5sum {output.cram:q} | awk '{{print $1}}' > {output.cram_md5:q}
        """


rule check_prevalidation:
    input:
        duplication_metrics=f"{OUTDIR}/{{sample}}/mark_duplicates/{{sample}}.metrics",
        agg_alignment_summary_metrics=f"{OUTDIR}/{{sample}}/collect_aggregation_metrics/{{sample}}.alignment_summary_metrics",
    output:
        duplication_value=f"{OUTDIR}/{{sample}}/check_prevalidation/{{sample}}.duplication_value.txt",
        chimerism_value=f"{OUTDIR}/{{sample}}/check_prevalidation/{{sample}}.chimerism_value.txt",
        is_outlier_data=f"{OUTDIR}/{{sample}}/check_prevalidation/{{sample}}.is_outlier_data.txt",
    log:
        f"{LOGDIR}/check_prevalidation/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/check_prevalidation/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/check_prevalidation.py \
            --duplication-metrics {input.duplication_metrics:q} \
            --alignment-summary-metrics {input.agg_alignment_summary_metrics:q} \
            --duplication-value {output.duplication_value:q} \
            --chimerism-value {output.chimerism_value:q} \
            --is-outlier-data {output.is_outlier_data:q}
        """


# MATE_NOT_FOUND is ignored since 2025-12-11, unlike WARP. EG100024 had exactly
# 150 reads whose mate was not found, and no other sample had any.
rule validate_sam_file:
    input:
        **REFERENCE,
        is_outlier_data=f"{OUTDIR}/{{sample}}/check_prevalidation/{{sample}}.is_outlier_data.txt",
        cram=f"{OUTDIR}/{{sample}}/convert_to_cram/{{sample}}.cram",
        crai=f"{OUTDIR}/{{sample}}/convert_to_cram/{{sample}}.cram.crai",
    output:
        f"{OUTDIR}/{{sample}}/validate_sam_file/{{sample}}.validation_report",
    params:
        java=java_options(),
    log:
        f"{LOGDIR}/validate_sam_file/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/validate_sam_file/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        picard \
            {params.java:q} \
            ValidateSamFile \
            -INPUT {input.cram:q} \
            -OUTPUT {output:q} \
            -REFERENCE_SEQUENCE {input.fasta:q} \
            -MAX_OUTPUT 1000000000 \
            -IGNORE MISSING_TAG_NM \
            -IGNORE MATE_NOT_FOUND \
            -MODE VERBOSE \
            -SKIP_MATE_VALIDATION "$(cat {input.is_outlier_data:q})" \
            -IS_BISULFITE_SEQUENCED false
        """
