rule validate_gvcf:
    input:
        gvcf = f"{config.results_dir}/reblock/{{sample}}.rb.g.vcf.gz",
        reference_fasta = config.input_files.reference_fasta,
        dbsnp_vcf = config.input_files.dbsnp_vcf,
        calling_interval_list = config.input_files.calling_interval_list,
    output:
        signal = f"{config.results_dir}/validate_gvcf/{{sample}}.OK",
    software: GATK_ENV
    log: f"{config.log_dir}/validate_gvcf/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        gatk --java-options "-Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m" \
        ValidateVariants \
            -V {input.gvcf} \
            -R {input.reference_fasta} \
            -L {input.calling_interval_list} \
            -gvcf \
            --validation-type-to-exclude ALLELES \
            --dbsnp {input.dbsnp_vcf} \
            --no-overlaps
        touch {output.signal}
        """
