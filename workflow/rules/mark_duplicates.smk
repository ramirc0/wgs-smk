import math

rule mark_duplicates:
    input:
        bam = f"{config.results_dir}/unmapped_bam_to_aligned/{{sample}}.bam"
    output:
        bam = f"{config.results_dir}/mark_duplicates/{{sample}}.bam",
        metrics = f"{config.results_dir}/mark_duplicates/{{sample}}.metrics"
    params:
        picard_tmp_dir = config.get("picard_tmp_dir", ".picard"),
        compression_level = config.get("compression_level", 2),
        mem_mb = lambda wildcards, resources: int(math.ceil(resources.mem_mb * 0.80))
    software: SOFTWARE_ENV
    log: f"{config.log_dir}/mark_duplicates/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1

        PICARD_TMP=$(mktemp -d -p {params.picard_tmp_dir} picard_tmp.XXXXXX)
        trap "rm -rf $PICARD_TMP" EXIT

        picard -Dsamjdk.compression_level={params.compression_level} \
        -Dorg.xerial.snappy.tempdir=$PICARD_TMP -Djava.io.tmpdir=$PICARD_TMP \
        -Xms{params.mem_mb}m -Xmx{params.mem_mb}m \
        MarkDuplicates \
            -INPUT {input.bam} \
            -OUTPUT {output.bam} \
            -METRICS_FILE {output.metrics} \
            -VALIDATION_STRINGENCY SILENT \
            -OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            -ASSUME_SORT_ORDER "queryname" \
            -CLEAR_DT "false" \
            -ADD_PG_TAG_TO_READS false
        """
