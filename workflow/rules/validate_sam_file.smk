rule validate_sam_file:
    input:
        duplication_value = f"{config.results_dir}/check_prevalidation/{{sample}}.duplication_value.txt",
        chimerism_value = f"{config.results_dir}/check_prevalidation/{{sample}}.chimerism_value.txt",
        bam = f"{config.results_dir}/convert_to_cram/{{sample}}.cram",
        reference_fasta = config.input_files.reference_fasta,
    output:
        validation_report = f"{config.results_dir}/validate_sam_file/{{sample}}.validation_report",
    software: SOFTWARE_ENV
    log: f"{config.log_dir}/validate_sam_file/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1

        CHIMERISM_THRESHOLD=0.15
        DUPLICATION_THRESHOLD=0.30 

        OBSERVED_DUPLICATION_VALUE=$(cat {input.duplication_value})
        OBSERVED_CHIMERISM_VALUE=$(cat {input.chimerism_value})

        if (( $(echo "$OBSERVED_DUPLICATION_VALUE > $DUPLICATION_THRESHOLD" | bc -l) )) || \
   (( $(echo "$OBSERVED_CHIMERISM_VALUE > $CHIMERISM_THRESHOLD" | bc -l) )); then
            echo "Skipping mate validation..."
            FLAG="-SKIP_MATE_VALIDATION true"
        else
            echo "Running mate validation..."
            FLAG="-SKIP_MATE_VALIDATION false"
        fi
        
        picard -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m \
        ValidateSamFile \
            -INPUT {input.bam} \
            -OUTPUT {output.validation_report} \
            -REFERENCE_SEQUENCE {input.reference_fasta} \
            -MAX_OUTPUT 1000000000 \
            -IGNORE MISSING_TAG_NM \
            -IGNORE MATE_NOT_FOUND \
            -MODE VERBOSE \
            $FLAG \
            -IS_BISULFITE_SEQUENCED false
        """

# Added MATE_NOT_FOUND to the ignore list (2025-12-11). EG100024 produced this error. Exactly 150 reads were marked as mate not found. This could be handled upstream
# at convert_to_cram by adding a step which includes "samtools fixmate -r", however given that EG100024 was the only sample to produce this error, it is 
# more than likely an edge case upstream of BSQR. Further testing is required to confirm this.
