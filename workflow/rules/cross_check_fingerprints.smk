rule cross_check_fingerprints:
    input:
        bam = f"{config.results_dir}/sort_bam/{{sample}}.bam",
        haplotype_database_file = config.input_files.haplotype_database_file
    output:
        metrics = f"{config.results_dir}/cross_check_fingerprints/{{sample}}.metrics"
    params:
        lod_threshold = config.get("lod_threshold", -20.0),
        cross_check_fingerprints_by = config.get("cross_check_fingerprints_by", "READGROUP")
    software: SOFTWARE_ENV
    log: f"{config.log_dir}/cross_check_fingerprints/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        picard -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 \
        -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m \
        CrosscheckFingerprints \
            -OUTPUT {output.metrics} \
            -HAPLOTYPE_MAP {input.haplotype_database_file} \
            -EXPECT_ALL_GROUPS_TO_MATCH true \
            -INPUT {input.bam} \
            -LOD_THRESHOLD {params.lod_threshold} \
            -CROSSCHECK_BY {params.cross_check_fingerprints_by}
        """
