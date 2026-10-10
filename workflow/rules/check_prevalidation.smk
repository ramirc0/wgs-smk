# TODO: Add some kind of check to return a boolean if duplication or chimerism values exceed these thresholds:
# Float max_duplication_in_reasonable_sample = 0.30
# Float max_chimerism_in_reasonable_sample = 0.15

rule check_prevalidation:
    input:
        duplication_metrics = f"{config.results_dir}/mark_duplicates/{{sample}}.metrics",
        agg_alignment_summary_metrics = f"{config.results_dir}/collect_aggregation_metrics/{{sample}}.alignment_summary_metrics",
    output:
        duplication_value = f"{config.results_dir}/check_prevalidation/{{sample}}.duplication_value.txt",
        chimerism_value = f"{config.results_dir}/check_prevalidation/{{sample}}.chimerism_value.txt"
    software: SOFTWARE_ENV
    log: f"{config.log_dir}/check_prevalidation/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        python {workflow.basedir}/scripts/check_prevalidation.py \
            --duplication-metrics {input.duplication_metrics} \
            --alignment-summary-metrics {input.agg_alignment_summary_metrics} \
            --duplication-value {output.duplication_value} \
            --chimerism-value {output.chimerism_value}
        """
