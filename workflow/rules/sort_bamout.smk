rule sort_bamout:
    input:
        bam = f"{config.results_dir}/haplotype_caller/{{sample}}/{{scatter}}.bamout.bam",
    output:
        bam = f"{config.results_dir}/sort_bamout/{{sample}}/{{scatter}}.bam",
        md5 = f"{config.results_dir}/sort_bamout/{{sample}}/{{scatter}}.bam.md5",
    software: SOFTWARE_ENV
    log: f"{config.log_dir}/sort_bamout/{{sample}}/{{scatter}}.log"
    shell:
        """
        exec >> {log} 2>&1
        samtools sort -@ {resources.threads} -m 2G -o {output.bam} {input.bam} && \
        samtools index {output.bam} && \
        md5sum {output.bam} > {output.md5}
        """
