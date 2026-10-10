def reformat_sequence_group_interval(wildcards):
    reformatted_intervals = ""
    temp_path = os.path.join(config.tmp_dir, f"{uuid.uuid4()}.tsv")
    shutil.copy2(SEQUENCE_GROUPS_PATH, temp_path)
    try:
        SEQUENCE_GROUPS = []
        with open(temp_path, "r") as tsv_file:
            all_content = tsv_file.read()
            for line in all_content.strip().split('\n'):
                line = line.strip()
                if line:
                    sequences = line.split("\t")
                    SEQUENCE_GROUPS.append(sequences)
    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)
    sequences = SEQUENCE_GROUPS[int(wildcards.group)]
    for sequence in sequences:
        reformatted_intervals += f"-L {sequence} "
    return reformatted_intervals

def get_bin_base_qualities_flags(_):
    reformatted_flags = ""
    if config.bin_base_qualities:
        for i in range(10, 40, 10):
            reformatted_flags += f"--static-quantized-quals {i} "
    return reformatted_flags

def get_bin_somatic_base_qualities_flags(_):
    reformatted_flags = ""
    if config.bin_base_qualities & config.somatic:
        for i in range(30, 60, 10):
            reformatted_flags += f"--static-quantized-quals {i} "
    return reformatted_flags

rule apply_bqsr:
    input:
        bam = f"{config.results_dir}/sort_bam/{{sample}}.bam",
        bam_index = f"{config.results_dir}/sort_bam/{{sample}}.bam.bai",
        reference_fasta = config.input_files.reference_fasta,
        bqsr_report = f"{config.results_dir}/gather_bqsr_reports/{{sample}}.recal_data.csv"
    output:
        bam = f"{config.results_dir}/apply_bqsr/{{sample}}/{{group}}.bam"
    params:
        compression_level = config.get("compression_level", 2),
        bin_base_qualities_flags = get_bin_base_qualities_flags,
        bin_somatic_base_qualities_flags = get_bin_somatic_base_qualities_flags,
        sequence_group_interval = reformat_sequence_group_interval
    software: GATK_ENV
    log: f"{config.log_dir}/apply_bqsr/{{sample}}/{{group}}.log"
    shell:
        """
        exec >> {log} 2>&1
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
        -Xlog:gc*:file={log}:time,uptime,level,tags \
        -Dsamjdk.compression_level={params.compression_level} \
        -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m" \
        ApplyBQSR \
            --create-output-bam-md5 \
            --add-output-sam-program-record \
            -R {input.reference_fasta} \
            -I {input.bam} \
            --use-original-qualities \
            -O {output.bam} \
            -bqsr {input.bqsr_report} \
            {params.bin_base_qualities_flags} \
            {params.bin_somatic_base_qualities_flags} \
            {params.sequence_group_interval}
        """
