rule convert_to_cram:
    input:
        bam = f"{config.results_dir}/gather_recalibrated_bams/{{sample}}.bam",
        reference_fasta = config.input_files.reference_fasta,
    output:
        cram_md5 = f"{config.results_dir}/convert_to_cram/{{sample}}.cram.md5",
        cram = f"{config.results_dir}/convert_to_cram/{{sample}}.cram"
    software: SOFTWARE_ENV
    log: f"{config.log_dir}/convert_to_cram/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1

        REF_CACHE_TMP=$(mktemp -d -p {config.tmp_dir} ref_cache.XXXXXX)
        trap "rm -rf $REF_CACHE_TMP" EXIT

        echo "Converting BAM to CRAM..."
        samtools view -@ {resources.threads} --output-fmt-option version=3.0 -C -T {input.reference_fasta} {input.bam} | \
        tee {output.cram} | \
        md5sum | awk '{{print $1}}' > {output.cram_md5}
        
        echo "Creating a REF_CACHE directory for CRAM indexing..."
        seq_cache_populate.pl -root $REF_CACHE_TMP {input.reference_fasta}
        export REF_PATH=:
        export REF_CACHE=$REF_CACHE_TMP/%2s/%2s/%s

        echo "Indexing CRAM..."
        samtools index -@ {resources.threads} {output.cram}
        """

