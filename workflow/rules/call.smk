# DRAGEN-mode variant calling to a reblocked, validated gVCF.


rule calibrate_dragen_str_model:
    input:
        **REFERENCE,
        bam=f"{OUTDIR}/{{sample}}/gather_recalibrated_bams/{{sample}}.bam",
        str_table_file=REFS["str_table_file"],
    output:
        f"{OUTDIR}/{{sample}}/calibrate_dragen_str_model/{{sample}}.dragstr",
    params:
        java=java_options(),
    log:
        f"{LOGDIR}/calibrate_dragen_str_model/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/calibrate_dragen_str_model/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        gatk --java-options \
        "{params.java} -Dgatk_stacktrace_on_user_exception=true -Dsamjdk.reference_fasta={input.fasta}" \
            CalibrateDragstrModel \
            -R {input.fasta:q} \
            -I {input.bam:q} \
            -str {input.str_table_file:q} \
            -O {output:q} \
            --threads {threads:q}
        """


rule haplotype_caller:
    input:
        **REFERENCE,
        bam=f"{OUTDIR}/{{sample}}/gather_recalibrated_bams/{{sample}}.bam",
        interval_list=lambda w: SCATTERED_INTERVAL_LIST.format(scatter=int(w.scatter)),
        contamination=f"{OUTDIR}/{{sample}}/check_contamination/{{sample}}.preBqsr.selfSM",
        dragen_str_model=f"{OUTDIR}/{{sample}}/calibrate_dragen_str_model/{{sample}}.dragstr",
    output:
        gvcf=f"{OUTDIR}/{{sample}}/haplotype_caller/{{scatter}}.g.vcf.gz",
        bam=f"{OUTDIR}/{{sample}}/haplotype_caller/{{scatter}}.bamout.bam",
    params:
        compression_level=lambda w: config.get("compression_level", 2),
        java=java_options(),
    log:
        f"{LOGDIR}/haplotype_caller/{{sample}}.{{scatter}}.txt",
    benchmark:
        f"{BENCHDIR}/haplotype_caller/{{sample}}.{{scatter}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        CONTAMINATION=$(awk 'NR==2 {{print $7}}' {input.contamination:q})

        gatk --java-options "{params.java} -Dsamjdk.compression_level={params.compression_level} -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            HaplotypeCaller \
            -R {input.fasta:q} \
            -I {input.bam:q} \
            -L {input.interval_list:q} \
            -O {output.gvcf:q} \
            -contamination "$CONTAMINATION" \
            -G StandardAnnotation \
            -G StandardHCAnnotation \
            -G AS_StandardAnnotation \
            --dragen-mode \
            --dragstr-params-path {input.dragen_str_model:q} \
            --native-pair-hmm-threads {threads:q} \
            -GQB 10 \
            -GQB 20 \
            -GQB 30 \
            -GQB 40 \
            -GQB 50 \
            -GQB 60 \
            -GQB 70 \
            -GQB 80 \
            -GQB 90 \
            -ERC GVCF \
            -bamout {output.bam:q}
        """


rule dragen_hardfilter_gvcf:
    input:
        f"{OUTDIR}/{{sample}}/haplotype_caller/{{scatter}}.g.vcf.gz",
    output:
        f"{OUTDIR}/{{sample}}/dragen_hardfilter_gvcf/{{scatter}}.hardfiltered.g.vcf.gz",
    params:
        java=java_options(),
    log:
        f"{LOGDIR}/dragen_hardfilter_gvcf/{{sample}}.{{scatter}}.txt",
    benchmark:
        f"{BENCHDIR}/dragen_hardfilter_gvcf/{{sample}}.{{scatter}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        gatk --java-options "{params.java}" \
            VariantFiltration \
            -V {input:q} \
            --filter-expression "QUAL < 10.4139" \
            --filter-name "DRAGENHardQUAL" \
            -O {output:q}
        """


rule merge_gvcfs:
    input:
        expand(
            f"{OUTDIR}/{{sample}}/dragen_hardfilter_gvcf/{{scatter}}.hardfiltered.g.vcf.gz",
            scatter=SCATTERS,
            allow_missing=True,
        ),
    output:
        gvcf=f"{OUTDIR}/{{sample}}/merge_gvcfs/{{sample}}.merged.g.vcf.gz",
        tbi=f"{OUTDIR}/{{sample}}/merge_gvcfs/{{sample}}.merged.g.vcf.gz.tbi",
    params:
        gvcfs=lambda w, input: repeat_flag("-INPUT", input),
        java=java_options(),
    log:
        f"{LOGDIR}/merge_gvcfs/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/merge_gvcfs/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        # The scatters are in genomic order, so GatherVcfs block-copies them.
        # Block copy cannot index, hence tabix.
        picard \
            {params.java:q} \
            GatherVcfs \
            {params.gvcfs:q} \
            -CREATE_INDEX false \
            -OUTPUT {output.gvcf:q}
        tabix -p vcf {output.gvcf:q}
        """


rule sort_bamout:
    input:
        f"{OUTDIR}/{{sample}}/haplotype_caller/{{scatter}}.bamout.bam",
    output:
        bam=f"{OUTDIR}/{{sample}}/sort_bamout/{{scatter}}.bam",
        md5=f"{OUTDIR}/{{sample}}/sort_bamout/{{scatter}}.bam.md5",
    log:
        f"{LOGDIR}/sort_bamout/{{sample}}.{{scatter}}.txt",
    benchmark:
        f"{BENCHDIR}/sort_bamout/{{sample}}.{{scatter}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        samtools sort -@ {threads:q} -m 2G -o {output.bam:q} {input:q}
        samtools index {output.bam:q}
        md5sum {output.bam:q} > {output.md5:q}
        """


rule merge_bamout:
    input:
        expand(
            f"{OUTDIR}/{{sample}}/sort_bamout/{{scatter}}.bam",
            scatter=SCATTERS,
            allow_missing=True,
        ),
    output:
        bam=f"{OUTDIR}/{{sample}}/merge_bamout/{{sample}}.bam",
        bam_index=f"{OUTDIR}/{{sample}}/merge_bamout/{{sample}}.bam.bai",
    log:
        f"{LOGDIR}/merge_bamout/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/merge_bamout/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        samtools merge -@ {threads:q} -o {output.bam:q} {input:q}
        samtools index -@ {threads:q} {output.bam:q}
        """


rule reblock:
    input:
        **REFERENCE,
        gvcf=f"{OUTDIR}/{{sample}}/merge_gvcfs/{{sample}}.merged.g.vcf.gz",
        tbi=f"{OUTDIR}/{{sample}}/merge_gvcfs/{{sample}}.merged.g.vcf.gz.tbi",
    output:
        f"{OUTDIR}/{{sample}}/reblock/{{sample}}.rb.g.vcf.gz",
    params:
        java=java_options(),
    log:
        f"{LOGDIR}/reblock/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/reblock/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        gatk --java-options "{params.java}" \
            ReblockGVCF \
            -R {input.fasta:q} \
            -V {input.gvcf:q} \
            -do-qual-approx \
            --floor-blocks \
            -GQB 20 \
            -GQB 30 \
            -GQB 40 \
            --add-site-filters-to-genotype \
            -O {output:q}
        """


rule validate_gvcf:
    input:
        **REFERENCE,
        gvcf=f"{OUTDIR}/{{sample}}/reblock/{{sample}}.rb.g.vcf.gz",
        dbsnp_vcf=REFS["dbsnp_vcf"],
        calling_interval_list=REFS["calling_interval_list"],
    output:
        f"{OUTDIR}/{{sample}}/validate_gvcf/{{sample}}.OK",
    params:
        java=java_options(),
    log:
        f"{LOGDIR}/validate_gvcf/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/validate_gvcf/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        gatk --java-options "{params.java}" \
            ValidateVariants \
            -V {input.gvcf:q} \
            -R {input.fasta:q} \
            -L {input.calling_interval_list:q} \
            -gvcf \
            --validation-type-to-exclude ALLELES \
            --dbsnp {input.dbsnp_vcf:q} \
            --no-overlaps
        touch {output:q}
        """
