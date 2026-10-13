# Reference indexes and the files built from them: the DRAGEN hash table, the
# STR table, reheadered interval lists and dbSNP, and the interval scatter.


# Writes beside the FASTA. ancient() keeps existing index files even when the
# FASTA is newer, as the parse-time indexing before this rule did.
rule reference_index:
    input:
        ancient(REFERENCE_FASTA),
    output:
        fasta_index=REFERENCE["fasta_index"],
        fasta_dict=REFERENCE["fasta_dict"],
    log:
        f"{LOGDIR}/reference_index.txt",
    benchmark:
        f"{BENCHDIR}/reference_index.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        samtools faidx {input:q}
        samtools dict {input:q} -o {output.fasta_dict:q}
        """


# Run-independent: built once under results/refs/, reused across runs.
rule dragen_index:
    input:
        **REFERENCE,
    output:
        directory(DRAGEN_REFERENCE),
    log:
        f"{LOGDIR}/dragen_index.txt",
    benchmark:
        f"{BENCHDIR}/dragen_index.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        mkdir -p {output:q}
        dragen-os \
            --build-hash-table true \
            --ht-reference {input.fasta:q} \
            --output-directory {output:q}
        """


rule compose_str_table:
    input:
        **REFERENCE,
    output:
        STR_TABLE,
    params:
        java=java_options(),
    log:
        f"{LOGDIR}/compose_str_table.txt",
    benchmark:
        f"{BENCHDIR}/compose_str_table.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        gatk --java-options "{params.java}" \
            ComposeSTRTableFile \
            -R {input.fasta:q} \
            -O {output:q}
        """


# The Broad VCF's header lists the contigs of Homo_sapiens_assembly38.
rule prepare_dbsnp:
    input:
        fasta_index=REFERENCE["fasta_index"],
        vcf=REFS["dbsnp_vcf"],
    output:
        vcf=DBSNP,
        tbi=f"{DBSNP}.tbi",
    log:
        f"{LOGDIR}/prepare_dbsnp.txt",
    benchmark:
        f"{BENCHDIR}/prepare_dbsnp.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bcftools reheader -f {input.fasta_index:q} -o {output.vcf:q} {input.vcf:q}
        tabix -p vcf {output.vcf:q}
        """


# Picard checks an interval list's header against the BAM's. The Broad lists
# carry the dictionary of Homo_sapiens_assembly38.
rule reheader_interval_list:
    input:
        fasta_dict=REFERENCE["fasta_dict"],
        interval_list=lambda w: REFS[f"{w.kind}_interval_list"],
    output:
        f"{REFDIR}/{{kind}}.interval_list",
    log:
        f"{LOGDIR}/reheader_interval_list/{{kind}}.txt",
    benchmark:
        f"{BENCHDIR}/reheader_interval_list/{{kind}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        {{ cat {input.fasta_dict:q}; grep -v '^@' {input.interval_list:q}; }} > {output:q}
        """


# Per run: the scatter count is a run parameter.
rule scatter_interval_list:
    input:
        REFS["calling_interval_list"],
    output:
        temp(expand(SCATTERED_INTERVAL_LIST, scatter=SCATTERS)),
    params:
        # Not a temp() directory output. With no consumer, Snakemake deletes it
        # and the lists inside at the next checkpoint.
        outdir=f"{OUTDIR}/scatter_interval_list",
        scatter_count=lambda w: config["scatter_interval_list"]["scatter_count"],
        break_bands_at_multiples_of=lambda w: config["scatter_interval_list"][
            "break_bands_at_multiples_of"
        ],
        java=java_options(),
    log:
        f"{LOGDIR}/scatter_interval_list.txt",
    benchmark:
        f"{BENCHDIR}/scatter_interval_list.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/scatter_interval_list.py \
            --java-options="{params.java}" \
            --scatter-count {params.scatter_count:q} \
            --break-bands-at-multiples-of {params.break_bands_at_multiples_of:q} \
            --interval-list {input:q} \
            --outdir {params.outdir:q}
        """
