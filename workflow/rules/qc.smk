# Read, alignment, contamination, fingerprint and variant calling QC.


rule collect_quality_yield_metrics:
    input:
        f"{OUTDIR}/{{sample}}/unmapped_bam/{{sample}}.bam",
    output:
        f"{OUTDIR}/{{sample}}/collect_quality_yield_metrics/{{sample}}.metrics",
    log:
        f"{LOGDIR}/collect_quality_yield_metrics/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/collect_quality_yield_metrics/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        picard \
            -Xms{resources.mem_mb:q}m \
            -Xmx{resources.mem_mb:q}m \
            CollectQualityYieldMetrics \
            -INPUT {input:q} \
            -OQ true \
            -OUTPUT {output:q}
        """


rule collect_unsorted_readgroup_bam_quality_metrics:
    input:
        f"{OUTDIR}/{{sample}}/unmapped_bam_to_aligned/{{sample}}.bam",
    output:
        insert_size_metrics=f"{OUTDIR}/{{sample}}/collect_unsorted_readgroup_bam_quality_metrics/{{sample}}.insert_size_metrics",
        insert_size_histogram=f"{OUTDIR}/{{sample}}/collect_unsorted_readgroup_bam_quality_metrics/{{sample}}.insert_size_histogram.pdf",
    params:
        prefix=f"{OUTDIR}/{{sample}}/collect_unsorted_readgroup_bam_quality_metrics/{{sample}}",
    log:
        f"{LOGDIR}/collect_unsorted_readgroup_bam_quality_metrics/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/collect_unsorted_readgroup_bam_quality_metrics/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        picard \
            -Xms{resources.mem_mb:q}m \
            -Xmx{resources.mem_mb:q}m \
            CollectMultipleMetrics \
            -INPUT {input:q} \
            -OUTPUT {params.prefix:q} \
            -ASSUME_SORTED true \
            -PROGRAM null \
            -PROGRAM CollectBaseDistributionByCycle \
            -PROGRAM CollectInsertSizeMetrics \
            -PROGRAM MeanQualityByCycle \
            -PROGRAM QualityScoreDistribution \
            -METRIC_ACCUMULATION_LEVEL null \
            -METRIC_ACCUMULATION_LEVEL ALL_READS
        """


rule cross_check_fingerprints:
    input:
        bam=f"{OUTDIR}/{{sample}}/sort_bam/{{sample}}.bam",
        haplotype_database_file=REFS["haplotype_database_file"],
    output:
        f"{OUTDIR}/{{sample}}/cross_check_fingerprints/{{sample}}.metrics",
    params:
        lod_threshold=lambda w: config.get("lod_threshold", -20.0),
        cross_check_fingerprints_by=lambda w: config.get(
            "cross_check_fingerprints_by", "READGROUP"
        ),
    log:
        f"{LOGDIR}/cross_check_fingerprints/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/cross_check_fingerprints/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        picard \
            -Dsamjdk.buffer_size=131072 \
            -XX:GCTimeLimit=50 \
            -XX:GCHeapFreeLimit=10 \
            -Xms{resources.mem_mb:q}m \
            -Xmx{resources.mem_mb:q}m \
            CrosscheckFingerprints \
            -OUTPUT {output:q} \
            -HAPLOTYPE_MAP {input.haplotype_database_file:q} \
            -EXPECT_ALL_GROUPS_TO_MATCH true \
            -INPUT {input.bam:q} \
            -LOD_THRESHOLD {params.lod_threshold:q} \
            -CROSSCHECK_BY {params.cross_check_fingerprints_by:q}
        """


rule check_contamination:
    input:
        **REFERENCE,
        bam=f"{OUTDIR}/{{sample}}/sort_bam/{{sample}}.bam",
        bam_index=f"{OUTDIR}/{{sample}}/sort_bam/{{sample}}.bam.bai",
        ud=REFS["contamination_sites"]["ud"],
        bed=REFS["contamination_sites"]["bed"],
        mu=REFS["contamination_sites"]["mu"],
    output:
        f"{OUTDIR}/{{sample}}/check_contamination/{{sample}}.preBqsr.selfSM",
    params:
        contamination_underestimation_factor=lambda w: config.get(
            "contamination_underestimation_factor", 0.75
        ),
    log:
        f"{LOGDIR}/check_contamination/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/check_contamination/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/check_contamination.py \
            --threads {threads:q} \
            --output {output:q} \
            --bam {input.bam:q} \
            --reference {input.fasta:q} \
            --ud {input.ud:q} \
            --bed {input.bed:q} \
            --mu {input.mu:q} \
            --underestimation-factor {params.contamination_underestimation_factor:q}
        """


rule collect_readgroup_bam_quality_metrics:
    input:
        **REFERENCE,
        bam=f"{OUTDIR}/{{sample}}/gather_recalibrated_bams/{{sample}}.bam",
    output:
        alignment_metrics=f"{OUTDIR}/{{sample}}/collect_readgroup_bam_quality_metrics/{{sample}}.alignment_summary_metrics",
        gc_bias_detail_metrics=f"{OUTDIR}/{{sample}}/collect_readgroup_bam_quality_metrics/{{sample}}.gc_bias.detail_metrics",
        gc_bias_summary_metrics=f"{OUTDIR}/{{sample}}/collect_readgroup_bam_quality_metrics/{{sample}}.gc_bias.summary_metrics",
        gc_bias_pdf=f"{OUTDIR}/{{sample}}/collect_readgroup_bam_quality_metrics/{{sample}}.gc_bias.pdf",
    params:
        prefix=f"{OUTDIR}/{{sample}}/collect_readgroup_bam_quality_metrics/{{sample}}",
    log:
        f"{LOGDIR}/collect_readgroup_bam_quality_metrics/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/collect_readgroup_bam_quality_metrics/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        picard \
            -Xms{resources.mem_mb:q}m \
            -Xmx{resources.mem_mb:q}m \
            CollectMultipleMetrics \
            -INPUT {input.bam:q} \
            -REFERENCE_SEQUENCE {input.fasta:q} \
            -OUTPUT {params.prefix:q} \
            -ASSUME_SORTED true \
            -PROGRAM null \
            -PROGRAM CollectAlignmentSummaryMetrics \
            -PROGRAM CollectGcBiasMetrics \
            -METRIC_ACCUMULATION_LEVEL null \
            -METRIC_ACCUMULATION_LEVEL READ_GROUP
        """


rule collect_aggregation_metrics:
    input:
        **REFERENCE,
        bam=f"{OUTDIR}/{{sample}}/gather_recalibrated_bams/{{sample}}.bam",
    output:
        alignment_summary_metrics=f"{OUTDIR}/{{sample}}/collect_aggregation_metrics/{{sample}}.alignment_summary_metrics",
        error_summary_metrics=f"{OUTDIR}/{{sample}}/collect_aggregation_metrics/{{sample}}.error_summary_metrics",
        bait_bias_detail_metrics=f"{OUTDIR}/{{sample}}/collect_aggregation_metrics/{{sample}}.bait_bias_detail_metrics",
        bait_bias_summary_metrics=f"{OUTDIR}/{{sample}}/collect_aggregation_metrics/{{sample}}.bait_bias_summary_metrics",
        gc_bias_detail_metrics=f"{OUTDIR}/{{sample}}/collect_aggregation_metrics/{{sample}}.gc_bias.detail_metrics",
        gc_bias_pdf=f"{OUTDIR}/{{sample}}/collect_aggregation_metrics/{{sample}}.gc_bias.pdf",
        gc_bias_summary_metrics=f"{OUTDIR}/{{sample}}/collect_aggregation_metrics/{{sample}}.gc_bias.summary_metrics",
        insert_size_histogram_pdf=f"{OUTDIR}/{{sample}}/collect_aggregation_metrics/{{sample}}.insert_size_histogram.pdf",
        insert_size_metrics=f"{OUTDIR}/{{sample}}/collect_aggregation_metrics/{{sample}}.insert_size_metrics",
        pre_adapter_detail_metrics=f"{OUTDIR}/{{sample}}/collect_aggregation_metrics/{{sample}}.pre_adapter_detail_metrics",
        pre_adapter_summary_metrics=f"{OUTDIR}/{{sample}}/collect_aggregation_metrics/{{sample}}.pre_adapter_summary_metrics",
        quality_distribution_pdf=f"{OUTDIR}/{{sample}}/collect_aggregation_metrics/{{sample}}.quality_distribution.pdf",
        quality_distribution_metrics=f"{OUTDIR}/{{sample}}/collect_aggregation_metrics/{{sample}}.quality_distribution_metrics",
    params:
        prefix=f"{OUTDIR}/{{sample}}/collect_aggregation_metrics/{{sample}}",
    log:
        f"{LOGDIR}/collect_aggregation_metrics/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/collect_aggregation_metrics/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        picard \
            -Xms{resources.mem_mb:q}m \
            -Xmx{resources.mem_mb:q}m \
            CollectMultipleMetrics \
            -INPUT {input.bam:q} \
            -REFERENCE_SEQUENCE {input.fasta:q} \
            -OUTPUT {params.prefix:q} \
            -ASSUME_SORTED true \
            -PROGRAM null \
            -PROGRAM CollectAlignmentSummaryMetrics \
            -PROGRAM CollectInsertSizeMetrics \
            -PROGRAM CollectSequencingArtifactMetrics \
            -PROGRAM QualityScoreDistribution \
            -PROGRAM CollectGcBiasMetrics \
            -METRIC_ACCUMULATION_LEVEL null \
            -METRIC_ACCUMULATION_LEVEL SAMPLE \
            -METRIC_ACCUMULATION_LEVEL LIBRARY
        """


# `INPUT= x` with a space is the legacy syntax, kept verbatim.
rule calculate_readgroup_checksum:
    input:
        f"{OUTDIR}/{{sample}}/gather_recalibrated_bams/{{sample}}.bam",
    output:
        f"{OUTDIR}/{{sample}}/calculate_readgroup_checksum/{{sample}}.bam.md5",
    log:
        f"{LOGDIR}/calculate_readgroup_checksum/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/calculate_readgroup_checksum/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        picard \
            -Xms{resources.mem_mb:q}m \
            -Xmx{resources.mem_mb:q}m \
            CalculateReadGroupChecksum \
            INPUT= {input:q} \
            OUTPUT= {output:q}
        """


rule collect_variant_calling_metrics:
    input:
        **REFERENCE,
        gvcf=f"{OUTDIR}/{{sample}}/reblock/{{sample}}.rb.g.vcf.gz",
        dbsnp_vcf=REFS["dbsnp_vcf"],
        chromsizes=CHROMSIZES,
        evaluation_interval_list=REFS["evaluation_interval_list"],
    output:
        variant_calling_detail_metrics=f"{OUTDIR}/{{sample}}/collect_variant_calling_metrics/{{sample}}.variant_calling_detail_metrics",
        variant_calling_summary_metrics=f"{OUTDIR}/{{sample}}/collect_variant_calling_metrics/{{sample}}.variant_calling_summary_metrics",
    params:
        prefix=f"{OUTDIR}/{{sample}}/collect_variant_calling_metrics/{{sample}}",
        tmp_dir=lambda w: config["tmp_dir"],
    log:
        f"{LOGDIR}/collect_variant_calling_metrics/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/collect_variant_calling_metrics/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        METRICS_TMP=$(mktemp -d -p {params.tmp_dir:q} variant_calling_metrics.XXXXXX)
        trap "rm -rf $METRICS_TMP" EXIT

        TMP_DBSNP_VCF_BGZ="$METRICS_TMP/dbsnp.vcf.gz"
        TMP_DBSNP_VCF_REHEADER="$METRICS_TMP/dbsnp.reheader.vcf.gz"

        echo "Fixing dbSNP VCF header..."
        bgzip -@ {threads:q} -c {input.dbsnp_vcf:q} > "$TMP_DBSNP_VCF_BGZ"
        bcftools index "$TMP_DBSNP_VCF_BGZ"
        bcftools reheader \
            -f {input.fasta_index:q} \
            -o "$TMP_DBSNP_VCF_REHEADER" \
            "$TMP_DBSNP_VCF_BGZ"
        bcftools index -t "$TMP_DBSNP_VCF_REHEADER"

        echo "Running CollectVariantCallingMetrics..."
        picard \
            -Xms{resources.mem_mb:q}m \
            -Xmx{resources.mem_mb:q}m \
            CollectVariantCallingMetrics \
            -INPUT {input.gvcf:q} \
            -OUTPUT {params.prefix:q} \
            --DBSNP "$TMP_DBSNP_VCF_REHEADER" \
            -SEQUENCE_DICTIONARY {input.fasta_dict:q} \
            -TARGET_INTERVALS {input.evaluation_interval_list:q} \
            -GVCF_INPUT true
        """
