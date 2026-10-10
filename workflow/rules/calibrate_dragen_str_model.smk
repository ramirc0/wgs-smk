# TDOO: Include a rule to makte the str_table_file from scratch (i.e. a difference reference genome)
# gatk --java-options "-Xmx8g" \
#     ComposeSTRTableFile \
#     -R reference.fasta \
#     -O my_custom.str

rule calibrate_dragen_str_model:
    input:
        bam = f"{config.results_dir}/gather_recalibrated_bams/{{sample}}.bam",
        reference_fasta = config.input_files.reference_fasta,
        str_table_file = config.input_files.str_table_file,
    output:
        dragen_str_model = f"{config.results_dir}/calibrate_dragen_str_model/{{sample}}.dragstr",
    params:
    software: GATK_ENV
    log: f"{config.log_dir}/calibrate_dragen_str_model/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        gatk --java-options \
        "-Xmx{resources.mem_mb}m -Dgatk_stacktrace_on_user_exception=true -Dsamjdk.reference_fasta={input.reference_fasta}" \
        CalibrateDragstrModel \
            -R {input.reference_fasta} \
            -I {input.bam} \
            -str {input.str_table_file} \
            -O {output.dragen_str_model} \
            --threads {resources.threads}
        """
