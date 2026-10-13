# FASTQs to per read group uBAMs, then a duplicate-marked, coordinate-sorted BAM.


# One uBAM per read group, `<flowcell>.<lane>`. A large one is also split into
# `<read group>/shard_<n>.bam`, aligned per shard and gathered.
rule unmapped_bam:
    input:
        r1=lambda w: R1_OF[w.sample],
        r2=lambda w: R2_OF[w.sample],
    output:
        temp(directory(f"{OUTDIR}/{{sample}}/unmapped_bam")),
    params:
        platform=lambda w: config.get("platform", "Generic"),
        compression_level=lambda w: config.get("compression_level", 2),
        cutoff_gib=lambda w: config["split_large_read_group"]["cutoff_gib"],
        reads_per_file=lambda w: config["split_large_read_group"]["reads_per_file"],
        # Per FastqToSam. One runs per thread.
        java=java_options(heap_mb=8000),
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
            --jobs {threads:q} \
            --sample {wildcards.sample:q} \
            --platform {params.platform:q} \
            --r1 {input.r1:q} \
            --r2 {input.r2:q} \
            --cutoff-gib {params.cutoff_gib:q} \
            --reads-per-file {params.reads_per_file:q} \
            --compression-level {params.compression_level:q} \
            --outdir {output:q}
        """


# The only checkpoint. It keeps the read group list after the uBAMs are gone.
checkpoint read_groups:
    input:
        f"{OUTDIR}/{{sample}}/unmapped_bam",
    output:
        f"{OUTDIR}/{{sample}}/read_groups/{{sample}}.tsv",
    log:
        f"{LOGDIR}/read_groups/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/read_groups/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        cp {input:q}/read_groups.tsv {output:q}
        """


rule unmapped_bam_to_aligned:
    input:
        **REFERENCE,
        dragen_reference=DRAGEN_REFERENCE,
        ubams=f"{OUTDIR}/{{sample}}/unmapped_bam",
    output:
        temp(f"{OUTDIR}/{{sample}}/unmapped_bam_to_aligned/{{unit}}.bam"),
    params:
        bam=lambda w, input: f"{input.ubams}/{w.unit}.bam",
        compression_level=lambda w: config.get("compression_level", 2),
        java=java_options(heap_mb=8000),
    log:
        f"{LOGDIR}/unmapped_bam_to_aligned/{{sample}}.{{unit}}.txt",
    benchmark:
        f"{BENCHDIR}/unmapped_bam_to_aligned/{{sample}}.{{unit}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        dragen-os \
            -b {params.bam:q} \
            -r {input.dragen_reference:q} \
            --num-threads {threads:q} \
            --interleaved=1 \
            --preserve-map-align-order true \
            | picard \
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
            -ALIGNED_BAM /dev/stdin \
            -UNMAPPED_BAM {params.bam:q} \
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
            -PROGRAM_GROUP_COMMAND_LINE "dragen-os -b {params.bam} -r {input.dragen_reference} --interleaved=1" \
            -PROGRAM_GROUP_NAME "dragen-os" \
            -UNMAPPED_READ_STRATEGY COPY_TO_TAG \
            -ALIGNER_PROPER_PAIR_FLAGS true \
            -UNMAP_CONTAMINANT_READS true \
            -ADD_PG_TAG_TO_READS false
        """


rule gather_unsorted_bam_files:
    input:
        lambda w: expand(
            f"{OUTDIR}/{w.sample}/unmapped_bam_to_aligned/{w.read_group}/{{shard}}.bam",
            shard=read_groups(w.sample)[w.read_group],
        ),
    output:
        temp(f"{OUTDIR}/{{sample}}/gather_unsorted_bam_files/{{read_group}}.bam"),
    params:
        bams=lambda w, input: repeat_flag("-INPUT", input),
        compression_level=lambda w: config.get("compression_level", 2),
        java=java_options(),
    log:
        f"{LOGDIR}/gather_unsorted_bam_files/{{sample}}.{{read_group}}.txt",
    benchmark:
        f"{BENCHDIR}/gather_unsorted_bam_files/{{sample}}.{{read_group}}.tsv",
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
            -CREATE_INDEX false \
            -CREATE_MD5_FILE false
        """


rule mark_duplicates:
    input:
        lambda w: aligned_read_groups(w.sample).values(),
    output:
        bam=temp(f"{OUTDIR}/{{sample}}/mark_duplicates/{{sample}}.bam"),
        metrics=f"{OUTDIR}/{{sample}}/mark_duplicates/{{sample}}.metrics",
    params:
        bams=lambda w, input: repeat_flag("-INPUT", input),
        compression_level=lambda w: config.get("compression_level", 2),
        # Tile, x and y are fields 5-7. Picard's default parser takes only 5- or
        # 7-field names, and MOHD's carry the UMI as an 8th.
        read_name_regex=lambda w: r"(?:[^:]*:){4}([0-9]+):([0-9]+):([0-9]+)(?::[^:]*)?",
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
            {params.bams:q} \
            -OUTPUT {output.bam:q} \
            -METRICS_FILE {output.metrics:q} \
            -VALIDATION_STRINGENCY SILENT \
            -READ_NAME_REGEX {params.read_name_regex:q} \
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
