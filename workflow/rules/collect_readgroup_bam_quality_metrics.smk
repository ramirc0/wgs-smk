rule collect_readgroup_bam_quality_metrics:
    input:
        bam = f"{config.results_dir}/gather_recalibrated_bams/{{sample}}.bam",
        reference_fasta = config.input_files.reference_fasta
    output:
        alignment_metrics = f"{config.results_dir}/collect_readgroup_bam_quality_metrics/{{sample}}.alignment_summary_metrics",
        gc_bias_detail_metrics = f"{config.results_dir}/collect_readgroup_bam_quality_metrics/{{sample}}.gc_bias.detail_metrics",
        gc_bias_summary_metrics = f"{config.results_dir}/collect_readgroup_bam_quality_metrics/{{sample}}.gc_bias.summary_metrics",
        gc_bias_pdf = f"{config.results_dir}/collect_readgroup_bam_quality_metrics/{{sample}}.gc_bias.pdf",
        signal = f"{config.results_dir}/collect_readgroup_bam_quality_metrics/{{sample}}.continue"
    params:
        output_prefix = f"{config.results_dir}/collect_readgroup_bam_quality_metrics/{{sample}}"
    software: SOFTWARE_ENV
    log: f"{config.log_dir}/collect_readgroup_bam_quality_metrics/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        picard -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m \
        CollectMultipleMetrics \
            -INPUT {input.bam} \
            -REFERENCE_SEQUENCE {input.reference_fasta} \
            -OUTPUT {params.output_prefix} \
            -ASSUME_SORTED true \
            -PROGRAM null \
            -PROGRAM CollectAlignmentSummaryMetrics \
            -PROGRAM CollectGcBiasMetrics \
            -METRIC_ACCUMULATION_LEVEL null \
            -METRIC_ACCUMULATION_LEVEL READ_GROUP
        touch {output.signal}
        """
