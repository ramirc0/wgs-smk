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
        java=java_options(),
    log:
        f"{LOGDIR}/base_recalibrator/{{sample}}.{{group}}.txt",
    benchmark:
        f"{BENCHDIR}/base_recalibrator/{{sample}}.{{group}}.tsv",
    software:
        GATK_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        gatk --java-options "{params.java} -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            BaseRecalibrator \
            -R {input.fasta:q} \
            -I {input.bam:q} \
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
        java=java_options(),
    log:
        f"{LOGDIR}/gather_bqsr_reports/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/gather_bqsr_reports/{{sample}}.tsv",
    software:
        GATK_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        gatk --java-options "{params.java}" \
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
        java=java_options(),
    log:
        f"{LOGDIR}/apply_bqsr/{{sample}}.{{group}}.txt",
    benchmark:
        f"{BENCHDIR}/apply_bqsr/{{sample}}.{{group}}.tsv",
    software:
        GATK_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        gatk --java-options "{params.java} -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 \
        -Dsamjdk.compression_level={params.compression_level}" \
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
        bams=lambda w, input: repeat_flag("-INPUT", input),
        java=java_options(),
    log:
        f"{LOGDIR}/gather_recalibrated_bams/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/gather_recalibrated_bams/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        picard \
            {params.java:q} \
            -Dsamjdk.compression_level={params.compression_level:q} \
            GatherBamFiles \
            {params.bams:q} \
            -OUTPUT {output:q} \
            -CREATE_INDEX true \
            -CREATE_MD5_FILE true
        """
