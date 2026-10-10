rule dragen_hardfilter_gvcf:
    input:
        gvcf = f"{config.results_dir}/haplotype_caller/{{sample}}/{{scatter}}.g.vcf.gz",
    output:
        gvcf = f"{config.results_dir}/dragen_hardfilter_gvcf/{{sample}}/{{scatter}}.hardfiltered.g.vcf.gz",
    software: GATK_ENV
    log: f"{config.log_dir}/dragen_hardfilter_gvcf/{{sample}}_{{scatter}}.log"
    shell:
        """
        exec >> {log} 2>&1
        gatk --java-options "-Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m" \
        VariantFiltration \
            -V {input.gvcf} \
            --filter-expression "QUAL < 10.4139" \
            --filter-name "DRAGENHardQUAL" \
            -O {output.gvcf}
        """
