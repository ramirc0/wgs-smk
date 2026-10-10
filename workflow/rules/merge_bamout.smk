def get_sorted_bams(wildcards):
    return expand(
        f"{config.results_dir}/sort_bamout/{wildcards.sample}/{{scatter}}.bam",
        scatter=range(1, config.scatter_interval_list.scatter_count + 1),
    )

def reformat_sorted_bams(wildcards):
    reformatted_sorted_bams = ""
    for scatter in range(1, config.scatter_interval_list.scatter_count + 1):
        reformatted_sorted_bams += f"{config.results_dir}/sort_bamout/{wildcards.sample}/{scatter}.bam "
    return reformatted_sorted_bams

rule merge_bamout:
    input:
        bams = get_sorted_bams
    output:
        bam = f"{config.results_dir}/merge_bamout/{{sample}}.bam",
        bam_index = f"{config.results_dir}/merge_bamout/{{sample}}.bam.bai",
    params:
        bams = reformat_sorted_bams
    software: SOFTWARE_ENV
    log: f"{config.log_dir}/merge_bamout/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        samtools merge -@ {resources.threads} -o {output.bam} {params.bams}
        samtools index -@ {resources.threads} {output.bam}
        """
