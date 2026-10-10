rule sort_bam:
    input:
        bam = f"{config.results_dir}/mark_duplicates/{{sample}}.bam"
    output:
        bam = f"{config.results_dir}/sort_bam/{{sample}}.bam",
        bam_index = f"{config.results_dir}/sort_bam/{{sample}}.bam.bai",
        md5 = f"{config.results_dir}/sort_bam/{{sample}}.bam.md5"
    software: SOFTWARE_ENV
    log: f"{config.log_dir}/sort_bam/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        samtools sort -@ {resources.threads} -m 2G -o {output.bam} {input.bam} && \
        samtools index {output.bam} && \
        md5sum {output.bam} > {output.md5}
        """
