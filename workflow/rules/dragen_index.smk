rule dragen_index:
    input:
        fasta = config.input_files.reference_fasta
    output:
        dragen_reference = directory(f"{config.results_dir}/dragen_reference"),
        signal = f"{config.results_dir}/dragen_reference/.continue"
    software: SOFTWARE_ENV
    log: f"{config.log_dir}/dragen_index.log"
    shell:
        """
        exec >> {log} 2>&1
        dragen-os --build-hash-table true --ht-reference {input.fasta} --output-directory {output.dragen_reference}
        touch {output.signal}
        """
