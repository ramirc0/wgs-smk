rule collect_variant_calling_metrics:
    input:
        gvcf = f"{config.results_dir}/reblock/{{sample}}.rb.g.vcf.gz",
        dbsnp_vcf = config.input_files.dbsnp_vcf,
        chromsizes = config.input_files.chromsizes,
        reference_fasta_index = config.input_files.reference_fasta_index,
        reference_fasta_dict = config.input_files.reference_fasta_dict,
        evaluation_interval_list = config.input_files.evaluation_interval_list,
    output:
        variant_calling_detail_metrics = f"{config.results_dir}/collect_variant_calling_metrics/{{sample}}.variant_calling_detail_metrics",
        variant_calling_summary_metrics = f"{config.results_dir}/collect_variant_calling_metrics/{{sample}}.variant_calling_summary_metrics",
    params:
        output_prefix = f"{config.results_dir}/collect_variant_calling_metrics/{{sample}}",
        tmp_dir = config.tmp_dir,
        compression_level = config.get("compression_level", 2),
    software: SOFTWARE_ENV
    log: f"{config.log_dir}/collect_variant_calling_metrics/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        
        METRICS_TMP=$(mktemp -d -p {params.tmp_dir} variant_calling_metrics.XXXXXX)
        trap "rm -rf $METRICS_TMP" EXIT

        TMP_DBSNP_VCF_BGZ="$METRICS_TMP/dbsnp.vcf.gz"
        TMP_DBSNP_VCF_REHEADER="$METRICS_TMP/dbsnp.reheader.vcf.gz"
        
        echo "Fixing dbSNP VCF header..."
        bgzip -@ {resources.threads} -c {input.dbsnp_vcf} > $TMP_DBSNP_VCF_BGZ
        bcftools index $TMP_DBSNP_VCF_BGZ
        bcftools reheader -f {input.reference_fasta_index} -o $TMP_DBSNP_VCF_REHEADER $TMP_DBSNP_VCF_BGZ
        bcftools index -t $TMP_DBSNP_VCF_REHEADER
        
        echo "Running CollectVariantCallingMetrics..."
        picard -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m \
        CollectVariantCallingMetrics \
            -INPUT {input.gvcf} \
            -OUTPUT {params.output_prefix} \
            -DBSNP $TMP_DBSNP_VCF_REHEADER \
            -SEQUENCE_DICTIONARY {input.reference_fasta_dict} \
            -TARGET_INTERVALS {input.evaluation_interval_list} \
            -GVCF_INPUT true
        """
