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

        python workflow/rules/scripts/unmapped_bam.py {resources.mem_mb} $PICARD_TMP {wildcards.sample} {params.platform} {input.R1} {input.R2} {output.bam}
        """
