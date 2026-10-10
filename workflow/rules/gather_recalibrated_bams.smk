def get_recalibrated_bams(wildcards):
    return expand(
        f"{config.results_dir}/apply_bqsr/{wildcards.sample}/{{group}}.bam",
        group=list(range(NUM_SEQUENCE_GROUPS))
    )

def reformat_recalibrated_bams(wildcards):
    reformatted_bams = ""
    bqsr_bams_template = f"{config.results_dir}/apply_bqsr/{wildcards.sample}/{{group}}.bam"
    for group in list(range(NUM_SEQUENCE_GROUPS)):
        reformatted_bams += f"-INPUT {bqsr_bams_template.format(group=group)} "
    return reformatted_bams

rule gather_recalibrated_bams:
    input:
        bqsr_bams = get_recalibrated_bams
    output:
        bam = f"{config.results_dir}/gather_recalibrated_bams/{{sample}}.bam"
    params:
        compression_level = config.get("compression_level", 2),
        picard_tmp_dir = config.get("picard_tmp_dir", ".picard"),
        bqsr_bams = reformat_recalibrated_bams
    software: SOFTWARE_ENV
    log: f"{config.log_dir}/gather_recalibrated_bams/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1

        PICARD_TMP=$(mktemp -d -p {params.picard_tmp_dir} picard_tmp.XXXXXX)
        trap "rm -rf $PICARD_TMP" EXIT

        picard -Dsamjdk.compression_level={params.compression_level} \
        -Dorg.xerial.snappy.tempdir=$PICARD_TMP -Djava.io.tmpdir=$PICARD_TMP \
        -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m \
        GatherBamFiles \
            {params.bqsr_bams} \
            -OUTPUT {output.bam} \
            -CREATE_INDEX true \
            -CREATE_MD5_FILE true
        """
