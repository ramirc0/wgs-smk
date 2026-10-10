rule reblock:
    input:
        reference_fasta = config.input_files.reference_fasta,
        gvcf = f"{config.results_dir}/merge_gvcfs/{{sample}}.merged.g.vcf.gz",
    output:
        gvcf = f"{config.results_dir}/reblock/{{sample}}.rb.g.vcf.gz",
    software: GATK_ENV
    log: f"{config.log_dir}/reblock/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        gatk --java-options "-Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m" \
        ReblockGVCF \
            -R {input.reference_fasta} \
            -V {input.gvcf} \
            -do-qual-approx \
            --floor-blocks -GQB 20 -GQB 30 -GQB 40 \
            --add-site-filters-to-genotype \
            -O {output.gvcf}
        """
