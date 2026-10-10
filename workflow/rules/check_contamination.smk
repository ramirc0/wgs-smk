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
        python ./workflow/rules/scripts/check_contamination.py {resources.threads} {output.metrics} {input.bam} {input.reference_fasta} {input.contamination_sites_ud} {input.contamination_sites_bed} {input.contamination_sites_mu} {params.contamination_underestimation_factor}
        """
