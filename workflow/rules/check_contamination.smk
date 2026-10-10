rule check_contamination:
    input:
        bam = f"{config.results_dir}/sort_bam/{{sample}}.bam",
        bam_index = f"{config.results_dir}/sort_bam/{{sample}}.bam.bai",
        reference_fasta = config.input_files.reference_fasta,
        contamination_sites_ud = config.contamination_sites.ud,
        contamination_sites_bed = config.contamination_sites.bed,
        contamination_sites_mu = config.contamination_sites.mu
    output:
        metrics = f"{config.results_dir}/check_contamination/{{sample}}.preBqsr.selfSM"
    params:
        contamination_underestimation_factor = config.get("contamination_underestimation_factor", 0.75),
    software: SOFTWARE_ENV
    log: f"{config.log_dir}/check_contamination/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        python {workflow.basedir}/scripts/check_contamination.py \
            --threads {resources.threads} \
            --output {output.metrics} \
            --bam {input.bam} \
            --reference {input.reference_fasta} \
            --ud {input.contamination_sites_ud} \
            --bed {input.contamination_sites_bed} \
            --mu {input.contamination_sites_mu} \
            --underestimation-factor {params.contamination_underestimation_factor}
        """
