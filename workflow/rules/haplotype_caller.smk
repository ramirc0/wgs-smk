def get_scatter_interval(wildcards):
    scatter = int(wildcards.scatter)
    return f"{config.results_dir}/scatter_interval_list/temp_{scatter:04d}_of_{config.scatter_interval_list.scatter_count}/{scatter}scattered.interval_list"

rule haplotype_caller:
    input:
        reference_fasta = config.input_files.reference_fasta,
        bam = f"{config.results_dir}/gather_recalibrated_bams/{{sample}}.bam",
        interval_list = get_scatter_interval,
        contamination = f"{config.results_dir}/check_contamination/{{sample}}.preBqsr.selfSM",
        dragen_str_model = f"{config.results_dir}/calibrate_dragen_str_model/{{sample}}.dragstr",
    output:
        gvcf = f"{config.results_dir}/haplotype_caller/{{sample}}/{{scatter}}.g.vcf.gz",
        bam = f"{config.results_dir}/haplotype_caller/{{sample}}/{{scatter}}.bamout.bam",
    params:
        picard_tmp_dir = config.get("picard_tmp_dir", ".picard"),
        compression_level = config.get("compression_level", 2),
    software: GATK_ENV
    log: f"{config.log_dir}/haplotype_caller/{{sample}}/{{scatter}}.log"
    shell:
        """
        exec >> {log} 2>&1

        PICARD_TMP=$(mktemp -d -p {params.picard_tmp_dir} picard_tmp.XXXXXX)
        trap "rm -rf $PICARD_TMP" EXIT
        
        CONTAMINATION=$(awk 'NR==2 {{print $7}}' {input.contamination})

        gatk --java-options "-Dorg.xerial.snappy.tempdir=$PICARD_TMP -Djava.io.tmpdir=$PICARD_TMP -Dsamjdk.compression_level={params.compression_level} -Xmx{resources.mem_mb}m -Xms{resources.mem_mb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        HaplotypeCaller \
            -R {input.reference_fasta} \
            -I {input.bam} \
            -L {input.interval_list} \
            -O {output.gvcf} \
            -contamination $CONTAMINATION \
            -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
            --dragen-mode \
            --dragstr-params-path {input.dragen_str_model} \
            --native-pair-hmm-threads {resources.threads} \
            -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
            -ERC GVCF \
            -bamout {output.bam} 
        """
