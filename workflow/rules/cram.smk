# CRAM conversion and its validation.


rule convert_to_cram:
    input:
        **REFERENCE,
        bam=f"{OUTDIR}/{{sample}}/gather_recalibrated_bams/{{sample}}.bam",
    output:
        cram_md5=f"{OUTDIR}/{{sample}}/convert_to_cram/{{sample}}.cram.md5",
        cram=f"{OUTDIR}/{{sample}}/convert_to_cram/{{sample}}.cram",
    params:
        tmp_dir=lambda w: config["tmp_dir"],
    log:
        f"{LOGDIR}/convert_to_cram/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/convert_to_cram/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        REF_CACHE_TMP=$(mktemp -d -p {params.tmp_dir:q} ref_cache.XXXXXX)
        trap "rm -rf $REF_CACHE_TMP" EXIT

        echo "Converting BAM to CRAM..."
        samtools view \
            -@ {threads:q} \
            --output-fmt-option version=3.0 \
            -C \
            -T {input.fasta:q} \
            {input.bam:q} \
            | tee {output.cram:q} \
            | md5sum \
            | awk '{{print $1}}' \
            > {output.cram_md5:q}

        echo "Creating a REF_CACHE directory for CRAM indexing..."
        seq_cache_populate.pl -root "$REF_CACHE_TMP" {input.fasta:q}
        export REF_PATH=:
        export REF_CACHE="$REF_CACHE_TMP/%2s/%2s/%s"

        echo "Indexing CRAM..."
        samtools index -@ {threads:q} {output.cram:q}
        """


rule check_prevalidation:
    input:
        duplication_metrics=f"{OUTDIR}/{{sample}}/mark_duplicates/{{sample}}.metrics",
        agg_alignment_summary_metrics=f"{OUTDIR}/{{sample}}/collect_aggregation_metrics/{{sample}}.alignment_summary_metrics",
    output:
        duplication_value=f"{OUTDIR}/{{sample}}/check_prevalidation/{{sample}}.duplication_value.txt",
        chimerism_value=f"{OUTDIR}/{{sample}}/check_prevalidation/{{sample}}.chimerism_value.txt",
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
            --chimerism-value {output.chimerism_value:q}
        """


# MATE_NOT_FOUND is ignored since 2025-12-11. EG100024 had exactly 150 reads
# whose mate was not found. `samtools fixmate -r` in convert_to_cram could fix
# it, but EG100024 is the only sample affected, so it is likely an edge case
# upstream of BQSR.
rule validate_sam_file:
    input:
        **REFERENCE,
        duplication_value=f"{OUTDIR}/{{sample}}/check_prevalidation/{{sample}}.duplication_value.txt",
        chimerism_value=f"{OUTDIR}/{{sample}}/check_prevalidation/{{sample}}.chimerism_value.txt",
        cram=f"{OUTDIR}/{{sample}}/convert_to_cram/{{sample}}.cram",
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

        CHIMERISM_THRESHOLD=0.15
        DUPLICATION_THRESHOLD=0.30

        OBSERVED_DUPLICATION_VALUE=$(cat {input.duplication_value:q})
        OBSERVED_CHIMERISM_VALUE=$(cat {input.chimerism_value:q})

        if (( $(echo "$OBSERVED_DUPLICATION_VALUE > $DUPLICATION_THRESHOLD" | bc -l) )) || \
            (( $(echo "$OBSERVED_CHIMERISM_VALUE > $CHIMERISM_THRESHOLD" | bc -l) )); then
            echo "Skipping mate validation..."
            FLAG="-SKIP_MATE_VALIDATION true"
        else
            echo "Running mate validation..."
            FLAG="-SKIP_MATE_VALIDATION false"
        fi

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
            $FLAG \
            -IS_BISULFITE_SEQUENCED false
        """
