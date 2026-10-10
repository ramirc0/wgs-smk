rule collect_aggregation_metrics:
    input:
        bam = f"{config.results_dir}/gather_recalibrated_bams/{{sample}}.bam",
        reference_fasta = config.input_files.reference_fasta,
    output:
        alignment_summary_metrics = f"{config['results_dir']}/collect_aggregation_metrics/{{sample}}.alignment_summary_metrics",
        error_summary_metrics = f"{config['results_dir']}/collect_aggregation_metrics/{{sample}}.error_summary_metrics", 
        bait_bias_detail_metrics = f"{config['results_dir']}/collect_aggregation_metrics/{{sample}}.bait_bias_detail_metrics",
        bait_bias_summary_metrics = f"{config['results_dir']}/collect_aggregation_metrics/{{sample}}.bait_bias_summary_metrics",
        gc_bias_detail_metrics = f"{config['results_dir']}/collect_aggregation_metrics/{{sample}}.gc_bias.detail_metrics",
        gc_bias_pdf = f"{config['results_dir']}/collect_aggregation_metrics/{{sample}}.gc_bias.pdf",
        gc_bias_summary_metrics = f"{config['results_dir']}/collect_aggregation_metrics/{{sample}}.gc_bias.summary_metrics",
        insert_size_histogram_pdf = f"{config['results_dir']}/collect_aggregation_metrics/{{sample}}.insert_size_histogram.pdf",
        insert_size_metrics = f"{config['results_dir']}/collect_aggregation_metrics/{{sample}}.insert_size_metrics",
        pre_adapter_detail_metrics = f"{config['results_dir']}/collect_aggregation_metrics/{{sample}}.pre_adapter_detail_metrics",
        pre_adapter_summary_metrics = f"{config['results_dir']}/collect_aggregation_metrics/{{sample}}.pre_adapter_summary_metrics",
        quality_distribution_pdf  = f"{config['results_dir']}/collect_aggregation_metrics/{{sample}}.quality_distribution.pdf",
        quality_distribution_metrics = f"{config['results_dir']}/collect_aggregation_metrics/{{sample}}.quality_distribution_metrics",
        signal = f"{config.results_dir}/collect_aggregation_metrics/{{sample}}.continue"
    params:
        output_prefix = f"{config.results_dir}/collect_aggregation_metrics/{{sample}}"
    software: SOFTWARE_ENV
    log: f"{config.log_dir}/collect_aggregation_metrics/{{sample}}.log"
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
            -PROGRAM CollectInsertSizeMetrics \
            -PROGRAM CollectSequencingArtifactMetrics \
            -PROGRAM QualityScoreDistribution \
            -PROGRAM CollectGcBiasMetrics \
            -METRIC_ACCUMULATION_LEVEL null \
            -METRIC_ACCUMULATION_LEVEL SAMPLE \
            -METRIC_ACCUMULATION_LEVEL LIBRARY
        touch {output.signal}
        """
