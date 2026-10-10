def get_bqsr_reports(wildcards):
    return expand(
        f"{config.results_dir}/base_recalibrator/{wildcards.sample}/{{group}}.recal_data.csv",
        group=list(range(NUM_SEQUENCE_GROUPS))
    )

def reformat_sequence_bqsr_reports(wildcards):
    reformatted_reports = ""
    bqsr_report_template = f"{config.results_dir}/base_recalibrator/{wildcards.sample}/{{group}}.recal_data.csv"
    for group in list(range(NUM_SEQUENCE_GROUPS)):
        reformatted_reports += f"-I {bqsr_report_template.format(group=group)} "
    return reformatted_reports

rule gather_bqsr_reports:
    input:
        bqsr_reports = get_bqsr_reports
    output:
        bqsr_report = f"{config.results_dir}/gather_bqsr_reports/{{sample}}.recal_data.csv"
    params:
        sequence_bqsr_reports = reformat_sequence_bqsr_reports
    software: GATK_ENV
    log: f"{config.log_dir}/gather_bqsr_reports/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        gatk --java-options "-Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m" \
            GatherBQSRReports \
            {params.sequence_bqsr_reports} \
            -O {output.bqsr_report}
          """
