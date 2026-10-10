
def reformat_known_indel_sites_vcfs(known_indel_sites_vcfs):
    reformatted_files = ""
    for known_indel_sites_vcf in known_indel_sites_vcfs:
        reformatted_files += f"--known-sites {known_indel_sites_vcf} "
    return reformatted_files

def reformat_sequence_group_interval(wildcards):
    reformatted_intervals = ""
    temp_path = os.path.join(config.tmp_dir, f"{uuid.uuid4()}.tsv")
    shutil.copy2(SEQUENCE_GROUPS_PATH, temp_path)

    try:
        SEQUENCE_GROUPS = []
        with open(temp_path, "r") as tsv_file:
            all_content = tsv_file.read()
            for line in all_content.strip().split('\n'):
                line = line.strip()
                if line:
                    sequences = line.split("\t")
                    SEQUENCE_GROUPS.append(sequences)
    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)

    # print("DEBUG")
    # print(f"Sequence groups: {SEQUENCE_GROUPS}")
    # print(f"Sequence groups path: {SEQUENCE_GROUPS_PATH}")
    # print(f"Number of sequence groups: {NUM_SEQUENCE_GROUPS}")

    sequences = SEQUENCE_GROUPS[int(wildcards.group)]
    for sequence in sequences:
        reformatted_intervals += f"-L {sequence} "
    return reformatted_intervals

# TODO: Make this more robust (add index files to list of inputs). This way the rule CANNOT be run if not provided by user.
# for now we hack away.

rule base_recalibration:
    input:
        bam = f"{config.results_dir}/sort_bam/{{sample}}.bam",
        bam_index = f"{config.results_dir}/sort_bam/{{sample}}.bam.bai",
        reference_fasta = config.input_files.reference_fasta,
        reference_fasta_index = config.input_files.reference_fasta_index,
        reference_fasta_dict = config.input_files.reference_fasta_dict,
        haplotype_database_file = config.input_files.haplotype_database_file,
        known_indel_sites_vcfs = config.input_files.known_indel_sites_vcfs,
        dbsnp_vcf = config.input_files.dbsnp_vcf,
    output:
        racalibration_report = f"{config.results_dir}/base_recalibrator/{{sample}}/{{group}}.recal_data.csv"
    params:
        known_indel_sites_vcfs_formatted = reformat_known_indel_sites_vcfs(config.input_files.known_indel_sites_vcfs),
        sequence_group_interval = reformat_sequence_group_interval,
        tmp_dir = config.tmp_dir
    software: GATK_ENV
    log: f"{config.log_dir}/base_recalibrator/{{sample}}/{{group}}.log"
    shell:
        """
        exec >> {log} 2>&1
        
        BQSR_TMP=$(mktemp -d -p {params.tmp_dir} bqsr_tmp.XXXXXX)
        trap "rm -rf $BQSR_TMP" EXIT

        cp {input.reference_fasta} \
        {input.reference_fasta_index} \
        {input.reference_fasta_dict} \
        {input.bam} \
        {input.bam_index} \
        $BQSR_TMP

        REFERNCE_FASTA=$BQSR_TMP/$(basename {input.reference_fasta})
        BAM=$BQSR_TMP/$(basename {input.bam})

        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
        -Xlog:gc*:file={log}:time,uptime,level,tags \
        -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m" \
        BaseRecalibrator \
            -R $REFERNCE_FASTA\
            -I $BAM \
            --use-original-qualities \
            -O {output.racalibration_report} \
            --known-sites {input.dbsnp_vcf} \
            {params.known_indel_sites_vcfs_formatted} \
            {params.sequence_group_interval}
        """
