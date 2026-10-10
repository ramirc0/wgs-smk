def get_fastq_files(wildcards):
      return {
        "R1": READ_LOOKUP[wildcards.sample][0],
        "R2": READ_LOOKUP[wildcards.sample][1]
    }

rule unmapped_bam:
    input: unpack(get_fastq_files)
    output:
        bam = f"{config.results_dir}/unmapped_bam/{{sample}}.bam"
    params:
        platform = config.get("platform", "Generic"),
        picard_tmp_dir = config.get("picard_tmp_dir", ".picard")
    software: SOFTWARE_ENV
    log: f"{config.log_dir}/unmapped_bam/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1

        PICARD_TMP=$(mktemp -d -p {params.picard_tmp_dir} picard_tmp.XXXXXX)
        trap "rm -rf $PICARD_TMP" EXIT

        python {workflow.basedir}/scripts/unmapped_bam.py \
            --mem-mb {resources.mem_mb} \
            --tmp-dir $PICARD_TMP \
            --sample {wildcards.sample} \
            --platform {params.platform} \
            --r1 {input.R1} \
            --r2 {input.R2} \
            --output {output.bam}
        """
