# FASTQs to a duplicate-marked, coordinate-sorted BAM.


rule unmapped_bam:
    input:
        r1=lambda w: R1_OF[w.sample],
        r2=lambda w: R2_OF[w.sample],
    output:
        f"{OUTDIR}/{{sample}}/unmapped_bam/{{sample}}.bam",
    params:
        platform=lambda w: config.get("platform", "Generic"),
        java=java_options(),
    log:
        f"{LOGDIR}/unmapped_bam/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/unmapped_bam/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/unmapped_bam.py \
            --java-options="{params.java}" \
            --tmp-dir "${{TMPDIR:-/tmp}}" \
            --sample {wildcards.sample:q} \
            --platform {params.platform:q} \
            --r1 {input.r1:q} \
            --r2 {input.r2:q} \
            --output {output:q}
        """


rule unmapped_bam_to_aligned:
    input:
        **REFERENCE,
        dragen_reference=DRAGEN_REFERENCE,
        bam=f"{OUTDIR}/{{sample}}/unmapped_bam/{{sample}}.bam",
    output:
        f"{OUTDIR}/{{sample}}/unmapped_bam_to_aligned/{{sample}}.bam",
    params:
        aligned_bam=lambda w: f"{config['tmp_dir']}/{w.sample}.aligned.bam",
        unmap_contaminant_reads=lambda w: config.get("unmap_contaminant_reads", False),
        compression_level=lambda w: config.get("compression_level", 2),
        java=java_options(heap_mb=8000),
    log:
        f"{LOGDIR}/unmapped_bam_to_aligned/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/unmapped_bam_to_aligned/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        dragen-os \
            -b {input.bam:q} \
            -r {input.dragen_reference:q} \
            --num-threads {threads:q} \
            --interleaved=1 \
            --preserve-map-align-order true \
            | samtools view -@ {threads:q} -h -O BAM - \
            > {params.aligned_bam:q}
        picard \
            {params.java:q} \
            -Dsamjdk.compression_level={params.compression_level:q} \
            MergeBamAlignment \
            -TMP_DIR "${{TMPDIR:-/tmp}}" \
            -VALIDATION_STRINGENCY SILENT \
            -EXPECTED_ORIENTATIONS FR \
            -ATTRIBUTES_TO_RETAIN X0 \
            -ATTRIBUTES_TO_REMOVE RG \
            -ATTRIBUTES_TO_REMOVE NM \
            -ATTRIBUTES_TO_REMOVE MD \
            -ALIGNED_BAM {params.aligned_bam:q} \
            -UNMAPPED_BAM {input.bam:q} \
            -OUTPUT {output:q} \
            -REFERENCE_SEQUENCE {input.fasta:q} \
            -PAIRED_RUN true \
            -SORT_ORDER "unsorted" \
            -IS_BISULFITE_SEQUENCE false \
            -ALIGNED_READS_ONLY false \
            -CLIP_ADAPTERS false \
            -MAX_RECORDS_IN_RAM 2000000 \
            -ADD_MATE_CIGAR true \
            -MAX_INSERTIONS_OR_DELETIONS -1 \
            -PRIMARY_ALIGNMENT_STRATEGY MostDistant \
            -PROGRAM_RECORD_ID "dragen-os" \
            -PROGRAM_GROUP_VERSION "$(dragen-os --version)" \
            -PROGRAM_GROUP_COMMAND_LINE "dragen-os -b {input.bam} -r {input.dragen_reference} --interleaved=1" \
            -PROGRAM_GROUP_NAME "dragen-os" \
            -UNMAPPED_READ_STRATEGY COPY_TO_TAG \
            -ALIGNER_PROPER_PAIR_FLAGS true \
            -UNMAP_CONTAMINANT_READS {params.unmap_contaminant_reads:q} \
            -ADD_PG_TAG_TO_READS false

        rm {params.aligned_bam:q}
        """


rule mark_duplicates:
    input:
        f"{OUTDIR}/{{sample}}/unmapped_bam_to_aligned/{{sample}}.bam",
    output:
        bam=f"{OUTDIR}/{{sample}}/mark_duplicates/{{sample}}.bam",
        metrics=f"{OUTDIR}/{{sample}}/mark_duplicates/{{sample}}.metrics",
    params:
        compression_level=lambda w: config.get("compression_level", 2),
        java=java_options(),
    log:
        f"{LOGDIR}/mark_duplicates/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/mark_duplicates/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        picard \
            {params.java:q} \
            -Dsamjdk.compression_level={params.compression_level:q} \
            MarkDuplicates \
            -TMP_DIR "${{TMPDIR:-/tmp}}" \
            -INPUT {input:q} \
            -OUTPUT {output.bam:q} \
            -METRICS_FILE {output.metrics:q} \
            -VALIDATION_STRINGENCY SILENT \
            -OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            -ASSUME_SORT_ORDER "queryname" \
            -CLEAR_DT "false" \
            -ADD_PG_TAG_TO_READS false
        """


rule sort_bam:
    input:
        f"{OUTDIR}/{{sample}}/mark_duplicates/{{sample}}.bam",
    output:
        bam=f"{OUTDIR}/{{sample}}/sort_bam/{{sample}}.bam",
        bam_index=f"{OUTDIR}/{{sample}}/sort_bam/{{sample}}.bam.bai",
        md5=f"{OUTDIR}/{{sample}}/sort_bam/{{sample}}.bam.md5",
    log:
        f"{LOGDIR}/sort_bam/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/sort_bam/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        samtools sort -@ {threads:q} -m 2G -o {output.bam:q} {input:q}
        samtools index {output.bam:q}
        md5sum {output.bam:q} > {output.md5:q}
        """
