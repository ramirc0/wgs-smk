rule calculate_readgroup_checksum:
    input:
        bam = f"{config.results_dir}/gather_recalibrated_bams/{{sample}}.bam",
    output:
        md5 = f"{config.results_dir}/calculate_readgroup_checksum/{{sample}}.bam.md5",
    software: SOFTWARE_ENV
    log: f"{config.log_dir}/calculate_readgroup_checksum/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        picard -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m \
        CalculateReadGroupChecksum \
            INPUT= {input.bam} \
            OUTPUT= {output.md5}
        """
