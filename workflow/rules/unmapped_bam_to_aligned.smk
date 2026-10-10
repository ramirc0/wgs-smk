rule unmapped_bam_to_aligned:
    input: 
        fasta = config.input_files.reference_fasta,
        dragen_reference = rules.dragen_index.output.dragen_reference,
        bam = f"{config.results_dir}/unmapped_bam/{{sample}}.bam"
    output:
        bam = f"{config.results_dir}/unmapped_bam_to_aligned/{{sample}}.bam"
    params:
        tmp_dir = config.tmp_dir,
        unmap_contaminant_reads = config.get("unmap_contaminant_reads", False),
        compression_level = config.get("compression_level", 2),
        picard_tmp_dir = config.get("picard_tmp_dir", ".picard")
    software: SOFTWARE_ENV
    log : f"{config.log_dir}/unmapped_bam_to_aligned/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1

        PICARD_TMP=$(mktemp -d -p {params.picard_tmp_dir} picard_tmp.XXXXXX)
        trap "rm -rf $PICARD_TMP" EXIT

        dragen-os -b {input.bam} -r {input.dragen_reference} --num-threads {resources.threads} --interleaved=1 --preserve-map-align-order true | samtools view -@ {resources.threads} -h -O BAM - > {params.tmp_dir}/{wildcards.sample}.aligned.bam
        picard -Dsamjdk.compression_level={params.compression_level} -Dorg.xerial.snappy.tempdir=$PICARD_TMP -Djava.io.tmpdir=$PICARD_TMP -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m MergeBamAlignment \
            -VALIDATION_STRINGENCY SILENT \
            -EXPECTED_ORIENTATIONS FR \
            -ATTRIBUTES_TO_RETAIN X0 \
            -ATTRIBUTES_TO_REMOVE RG \
            -ATTRIBUTES_TO_REMOVE NM \
            -ATTRIBUTES_TO_REMOVE MD \
            -ALIGNED_BAM {params.tmp_dir}/{wildcards.sample}.aligned.bam \
            -UNMAPPED_BAM {input.bam} \
            -OUTPUT {output.bam} \
            -REFERENCE_SEQUENCE {input.fasta} \
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
            -UNMAP_CONTAMINANT_READS {params.unmap_contaminant_reads} \
            -ADD_PG_TAG_TO_READS false

        rm {params.tmp_dir}/{wildcards.sample}.aligned.bam
        """
