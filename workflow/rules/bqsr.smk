# Base quality score recalibration, scattered over sequence groups.


rule base_recalibration:
    input:
        **REFERENCE,
        bam=f"{OUTDIR}/{{sample}}/sort_bam/{{sample}}.bam",
        bam_index=f"{OUTDIR}/{{sample}}/sort_bam/{{sample}}.bam.bai",
        haplotype_database_file=REFS["haplotype_database_file"],
        known_indel_sites_vcfs=REFS["known_indel_sites_vcfs"],
        dbsnp_vcf=REFS["dbsnp_vcf"],
    output:
        f"{OUTDIR}/{{sample}}/base_recalibrator/{{group}}.recal_data.csv",
    params:
        known_sites=lambda w, input: repeat_flag(
            "--known-sites", input.known_indel_sites_vcfs
        ),
        intervals=lambda w: repeat_flag("-L", SEQUENCE_GROUPS[int(w.group)]),
        tmp_dir=lambda w: config["tmp_dir"],
    log:
        f"{LOGDIR}/base_recalibrator/{{sample}}.{{group}}.txt",
    benchmark:
        f"{BENCHDIR}/base_recalibrator/{{sample}}.{{group}}.tsv",
    software:
        GATK_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        BQSR_TMP=$(mktemp -d -p {params.tmp_dir:q} bqsr_tmp.XXXXXX)
        trap "rm -rf $BQSR_TMP" EXIT

        cp \
            {input.fasta:q} \
            {input.fasta_index:q} \
            {input.fasta_dict:q} \
            {input.bam:q} \
            {input.bam_index:q} \
            "$BQSR_TMP"

        REFERENCE_FASTA="$BQSR_TMP/$(basename {input.fasta:q})"
        BAM="$BQSR_TMP/$(basename {input.bam:q})"

        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
        -Xlog:gc*:file={log}:time,uptime,level,tags \
        -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m" \
            BaseRecalibrator \
            -R "$REFERENCE_FASTA" \
            -I "$BAM" \
            --use-original-qualities \
            -O {output:q} \
            --known-sites {input.dbsnp_vcf:q} \
            {params.known_sites:q} \
            {params.intervals:q}
        """


rule gather_bqsr_reports:
    input:
        expand(
            f"{OUTDIR}/{{sample}}/base_recalibrator/{{group}}.recal_data.csv",
            group=GROUPS,
            allow_missing=True,
        ),
    output:
        f"{OUTDIR}/{{sample}}/gather_bqsr_reports/{{sample}}.recal_data.csv",
    params:
        reports=lambda w, input: repeat_flag("-I", input),
    log:
        f"{LOGDIR}/gather_bqsr_reports/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/gather_bqsr_reports/{{sample}}.tsv",
    software:
        GATK_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        gatk --java-options "-Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m" \
            GatherBQSRReports \
            {params.reports:q} \
            -O {output:q}
        """


rule apply_bqsr:
    input:
        **REFERENCE,
        bam=f"{OUTDIR}/{{sample}}/sort_bam/{{sample}}.bam",
        bam_index=f"{OUTDIR}/{{sample}}/sort_bam/{{sample}}.bam.bai",
        bqsr_report=f"{OUTDIR}/{{sample}}/gather_bqsr_reports/{{sample}}.recal_data.csv",
    output:
        f"{OUTDIR}/{{sample}}/apply_bqsr/{{group}}.bam",
    params:
        compression_level=lambda w: config.get("compression_level", 2),
        bin_base_qualities_flags=lambda w: (
            repeat_flag("--static-quantized-quals", [10, 20, 30])
            if config["bin_base_qualities"]
            else []
        ),
        bin_somatic_base_qualities_flags=lambda w: (
            repeat_flag("--static-quantized-quals", [30, 40, 50])
            if config["bin_base_qualities"] and config["somatic"]
            else []
        ),
        intervals=lambda w: repeat_flag("-L", SEQUENCE_GROUPS[int(w.group)]),
    log:
        f"{LOGDIR}/apply_bqsr/{{sample}}.{{group}}.txt",
    benchmark:
        f"{BENCHDIR}/apply_bqsr/{{sample}}.{{group}}.tsv",
    software:
        GATK_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
        -Xlog:gc*:file={log}:time,uptime,level,tags \
        -Dsamjdk.compression_level={params.compression_level} \
        -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m" \
            ApplyBQSR \
            --create-output-bam-md5 \
            --add-output-sam-program-record \
            -R {input.fasta:q} \
            -I {input.bam:q} \
            --use-original-qualities \
            -O {output:q} \
            -bqsr {input.bqsr_report:q} \
            {params.bin_base_qualities_flags:q} \
            {params.bin_somatic_base_qualities_flags:q} \
            {params.intervals:q}
        """


rule gather_recalibrated_bams:
    input:
        expand(
            f"{OUTDIR}/{{sample}}/apply_bqsr/{{group}}.bam",
            group=GROUPS,
            allow_missing=True,
        ),
    output:
        f"{OUTDIR}/{{sample}}/gather_recalibrated_bams/{{sample}}.bam",
    params:
        compression_level=lambda w: config.get("compression_level", 2),
        picard_tmp_dir=lambda w: config.get("picard_tmp_dir", ".picard"),
        bams=lambda w, input: repeat_flag("-INPUT", input),
    log:
        f"{LOGDIR}/gather_recalibrated_bams/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/gather_recalibrated_bams/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        PICARD_TMP=$(mktemp -d -p {params.picard_tmp_dir:q} picard_tmp.XXXXXX)
        trap "rm -rf $PICARD_TMP" EXIT

        picard \
            -Dsamjdk.compression_level={params.compression_level:q} \
            -Dorg.xerial.snappy.tempdir="$PICARD_TMP" \
            -Djava.io.tmpdir="$PICARD_TMP" \
            -Xms{resources.mem_mb:q}m \
            -Xmx{resources.mem_mb:q}m \
            GatherBamFiles \
            {params.bams:q} \
            -OUTPUT {output:q} \
            -CREATE_INDEX true \
            -CREATE_MD5_FILE true
        """
