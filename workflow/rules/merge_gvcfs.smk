def get_gvcfs(wildcards):
    return expand(
        f"{config.results_dir}/dragen_hardfilter_gvcf/{wildcards.sample}/{{scatter}}.hardfiltered.g.vcf.gz",
        scatter=range(1, config.scatter_interval_list.scatter_count + 1),
    )

def reformat_gvcfs(wildcards):
    reformatted_vcfs = ""
    for scatter in range(1, config.scatter_interval_list.scatter_count + 1):
        reformatted_vcfs += f"-INPUT {config.results_dir}/dragen_hardfilter_gvcf/{wildcards.sample}/{scatter}.hardfiltered.g.vcf.gz "
    return reformatted_vcfs

rule merge_gvcfs:
    input:
        gvcfs = get_gvcfs
    output:
        gvcf = f"{config.results_dir}/merge_gvcfs/{{sample}}.merged.g.vcf.gz"
    params:
        gvcfs = reformat_gvcfs
    software: SOFTWARE_ENV
    log: f"{config.log_dir}/merge_gvcfs/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        picard -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m \
        MergeVcfs \
            {params.gvcfs} \
            -OUTPUT {output.gvcf}
        """
