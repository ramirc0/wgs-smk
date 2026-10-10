# WGS Snakemake Pipeline

## Overview

Here we present a [snakemake](https://github.com/snakemake/snakemake.git) pipeline for whole genome sequencing (WGS) analysis of germline single samples. Starting from paired-end FASTQ files, this pipeline performs alignment, duplicate marking, base quality score recalibration (BQSR), variant calling, and produces analysis-ready gVCFs suitable for downstream joint genotyping.

> [!NOTE]
> This pipeline is based on the [WARP Whole Genome Germline Single Sample pipeline](https://broadinstitute.github.io/warp/docs/Pipelines/Whole_Genome_Germline_Single_Sample/README) by the Broad Institute, [originally licensed](https://github.com/broadinstitute/warp/blob/develop/LICENSE) under BSD. See [THIRD_PARTY_LICENSES](THIRD_PARTY_LICENSES) for details.

## Example Usage

### Prerequisites

This pipeline takes paired-end FASTQ files as input and produces analysis-ready gVCFs. The output gVCFs are suitable as input for the [joint-genotyping-smk](https://github.com/christian728504/joint-genotyping-smk.git) pipeline for joint genotyping across multiple samples.

### Dependencies

This pipeline requires the following dependencies:

- [`conda` (preferably `mamba`)](https://github.com/conda-forge/miniforge)
- [`singularity`](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)
- (Optional) [`slurm`](https://slurm.schedmd.com/quickstart.html)

To start, first setup the runtime environment:

```bash
mamba env create --name=wgs-smk --file=environments/wgs-smk.yaml
mamba activate wgs-smk
```

Then, download the required resources:

```bash
curl -fsSL https://users.wenglab.org/ramirezc/mohd/resources/wgs-smk.tar.zst | tar --use-compress-program=zstd -xf -
```

### Input Metadata

This pipeline expects a tab-delimited metadata file with a header row. The file must contain at minimum the following columns:

- `Sample`: Unique sample identifiers.
- `R1`: Path to the forward (R1) FASTQ file.
- `R2`: Path to the reverse (R2) FASTQ file.

An example metadata file:

```
Sample	R1	R2
EG100001	/path/to/EG100001_R1.fastq.gz	/path/to/EG100001_R2.fastq.gz
EG100002	/path/to/EG100002_R1.fastq.gz	/path/to/EG100002_R2.fastq.gz
```

### Example Configfile

For a complete example, see [config/config.yml](config/config.yml). From this file you'll want to pay attention to the following sections:

```yaml
input_files:
  metadata: "jobs/2026-01-16/metadata.tsv"
  reference_fasta: "jobs/2026-01-16/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
```

- `input_files.metadata`: Path to the tab-delimited metadata file containing sample IDs and FASTQ paths.
- `input_files.reference_fasta`: Path to the reference FASTA file. The pipeline will automatically generate the `.fai` index, `.dict` dictionary, and `.chromsizes` file if they do not already exist.

```yaml
platform: "Illumina"
compression_level: 2
unmap_contaminant_reads: True
```

- `platform`: Sequencing platform (used in read group metadata). Default is `Illumina`.
- `compression_level`: BAM compression level. Default is `2`.
- `unmap_contaminant_reads`: Whether to unmap contaminant reads during alignment merging. Default is `True`.

```yaml
contamination_sites:
  ud: "jobs/2026-01-16/1000g.phase3.100k.b38.vcf.gz.dat.UD"
  bed: "jobs/2026-01-16/1000g.phase3.100k.b38.vcf.gz.dat.bed"
  mu: "jobs/2026-01-16/1000g.phase3.100k.b38.vcf.gz.dat.mu"
contamination_underestimation_factor: 0.75
```

- `contamination_sites`: Paths to VerifyBamID resource files for contamination estimation.
- `contamination_underestimation_factor`: Factor applied to contamination estimates. Default is `0.75`.

```yaml
scatter_interval_list:
  scatter_count: 10
  break_bands_at_multiples_of: 100000
```

- `scatter_interval_list.scatter_count`: Number of intervals to scatter for parallelized variant calling with HaplotypeCaller. Higher values increase parallelism.
- `scatter_interval_list.break_bands_at_multiples_of`: Break bands at multiples of this value when scattering intervals. Default is `100000`.

```yaml
picard_tmp_dir: ".picard"
local_tmp_dir: ".tmp"
tmp_dir: "/tmp"
```

- `picard_tmp_dir`: Temporary directory for Picard tools. In some compute environments, `/tmp` is mounted with `noexec`, which will cause Picard to fail. In this case, set this to a local directory such as `.picard`.
- `local_tmp_dir`: Local temporary directory for intermediate files.
- `tmp_dir`: General temporary directory for intermediate files.

```yaml
environments:
  default: "docker://clarity001/wgs-smk:latest"
  gatk: "docker://us.gcr.io/broad-gatk/gatk:latest"
  picard: "docker://broadinstitute/picard:latest"
```

- `environments.default`: The default container image. This is a custom image containing samtools, bcftools, picard, DRAGEN-OS, VerifyBamID, and other tools.
- `environments.gatk`: The GATK container image, used for GATK-specific rules.
- `environments.picard`: The Picard container image.

### Defining Resources

All resource definitions for rules are located in [profile/slurm/config.yaml](profile/slurm/config.yaml).

```yaml
  haplotype_caller:
    runtime: 720m
    constraint: cascadelake
    slurm_extra: "--exclude=z[1024,1062,1068] --cpu-freq=High-High:Performance"
    slurm_partition: 12hours,5days
    threads: 26
    cpus_per_task: 26
    mem: 250000MB
```

In this instance, we are requesting a job on a SLURM cluster with a runtime of 720 minutes, a specific constraint (i.e. cascadelake), specific partitions (i.e. 12hours or 5days), 26 threads, 26 CPUs per task, and 250000MB of memory.

These resource definitions are specific to the Weng Lab's SLURM cluster. If you are running the pipeline on your own cluster, you'll need to adjust these values accordingly. We recommend consulting the [snakemake-executor-plugin-slurm documentation](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) for further information.

### Pipeline Steps

The pipeline executes the following steps for each sample:

1. **DRAGEN Index** — Build a DRAGEN hash table from the reference FASTA.
2. **Unmapped BAM** — Convert paired-end FASTQs to an unmapped BAM (uBAM) with read group metadata extracted from FASTQ headers.
3. **Quality Yield Metrics** — Collect quality yield metrics on the uBAM.
4. **Alignment** — Align reads to the reference genome using DRAGEN-OS and merge with the uBAM via Picard `MergeBamAlignment`.
5. **Unsorted Read Group BAM Quality Metrics** — Collect insert size and other quality metrics on the unsorted aligned BAM.
6. **Mark Duplicates** — Mark duplicate reads using Picard `MarkDuplicates`.
7. **Sort BAM** — Coordinate-sort and index the deduplicated BAM.
8. **Cross-Check Fingerprints** — Verify sample identity consistency across read groups.
9. **Check Contamination** — Estimate cross-sample contamination using VerifyBamID.
10. **Base Recalibration** — Run GATK `BaseRecalibrator` across sequence groups in parallel.
11. **Gather BQSR Reports** — Merge per-group recalibration reports.
12. **Apply BQSR** — Apply base quality score recalibration across sequence groups in parallel.
13. **Gather Recalibrated BAMs** — Merge per-group recalibrated BAMs into a single BAM.
14. **Read Group BAM Quality Metrics** — Collect alignment summary and GC bias metrics on the final BAM.
15. **Aggregation Metrics** — Collect comprehensive QC metrics (alignment, insert size, sequencing artifacts, GC bias, quality distribution).
16. **Read Group Checksum** — Calculate an MD5 checksum for the final BAM.
17. **Convert to CRAM** — Convert the final BAM to CRAM format for archival storage.
18. **Check Prevalidation** — Extract duplication and chimerism values to determine validation stringency.
19. **Validate SAM File** — Validate the CRAM file with Picard `ValidateSamFile`.
20. **Calibrate DRAGEN STR Model** — Calibrate the DRAGEN short tandem repeat model using GATK `CalibrateDragstrModel`.
21. **Scatter Interval List** — Scatter calling intervals for parallelized variant calling.
22. **HaplotypeCaller** — Call germline variants per scatter interval in DRAGEN mode, producing per-interval gVCFs and BAMout files.
23. **DRAGEN Hard Filter gVCF** — Apply DRAGEN hard quality filtering to per-interval gVCFs.
24. **Merge gVCFs** — Merge per-interval hard-filtered gVCFs into a single sample gVCF.
25. **Sort BAMout / Merge BAMout** — Sort and merge per-interval BAMout files.
26. **Reblock** — Reblock the merged gVCF using GATK `ReblockGVCF` to reduce file size and standardize GQ bands.
27. **Validate gVCF** — Validate the final reblocked gVCF.
28. **Collect Variant Calling Metrics** — Collect variant calling detail and summary metrics.

### Running the Pipeline

With the hard part out of the way, you can now run the pipeline:

```bash
snakemake --workflow-profile profile/slurm
```

> [!NOTE]
> It is recommended that you run the pipeline in dry run mode first (add the `-n` flag). Also note that `snakemake` must be ran in the root of this repository.

This will run the pipeline and produce all outputs to the directory specified by `results_dir` in the config file (default: `results`).

### Outputs

The pipeline produces results organized under the `results/` directory:

| Directory | Key Outputs |
|---|---|
| `results/dragen_reference/` | DRAGEN hash table index for the reference genome |
| `results/unmapped_bam/` | Unmapped BAMs (uBAMs) with read group metadata |
| `results/collect_quality_yield_metrics/` | Quality yield metrics for each uBAM |
| `results/unmapped_bam_to_aligned/` | Aligned BAMs (merged with uBAM) |
| `results/collect_unsorted_readgroup_bam_quality_metrics/` | Insert size metrics and histogram PDFs |
| `results/mark_duplicates/` | Deduplicated BAMs, duplication metrics |
| `results/sort_bam/` | Coordinate-sorted BAMs, BAM indexes, MD5 checksums |
| `results/cross_check_fingerprints/` | Fingerprint cross-check metrics |
| `results/check_contamination/` | VerifyBamID `.selfSM` contamination estimates |
| `results/base_recalibrator/` | Per-sequence-group BQSR recalibration tables |
| `results/gather_bqsr_reports/` | Merged BQSR recalibration reports |
| `results/apply_bqsr/` | Per-sequence-group recalibrated BAMs |
| `results/gather_recalibrated_bams/` | Final merged recalibrated BAMs |
| `results/collect_readgroup_bam_quality_metrics/` | Read group alignment summary and GC bias metrics |
| `results/collect_aggregation_metrics/` | Comprehensive QC: alignment, insert size, GC bias, sequencing artifacts, quality distribution |
| `results/calculate_readgroup_checksum/` | BAM MD5 checksums |
| `results/convert_to_cram/` | CRAM files and CRAM MD5 checksums |
| `results/check_prevalidation/` | Duplication and chimerism value files |
| `results/validate_sam_file/` | Picard `ValidateSamFile` validation reports |
| `results/calibrate_dragen_str_model/` | DRAGEN STR model `.dragstr` files |
| `results/scatter_interval_list/` | Scattered interval lists for parallel variant calling |
| `results/haplotype_caller/` | Per-interval gVCFs (`.g.vcf.gz`) and BAMout files |
| `results/dragen_hardfilter_gvcf/` | DRAGEN hard-filtered per-interval gVCFs |
| `results/merge_gvcfs/` | Merged per-sample gVCFs |
| `results/sort_bamout/` | Sorted per-interval BAMout files |
| `results/merge_bamout/` | Merged BAMout BAMs |
| `results/reblock/` | Reblocked gVCFs (`.rb.g.vcf.gz`), ready for joint genotyping |
| `results/validate_gvcf/` | gVCF validation signals |
| `results/collect_variant_calling_metrics/` | Variant calling detail and summary metrics |

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.

> https://github.com/broadinstitute/warp/blob/develop/pipelines/wdl/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl

> https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/UnmappedBamToAlignedBam.wdl

> https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/BamToGvcf.wdl

## Questions

If you have any questions or would like to provide constructive feedback, please open an issue or reach out to the MOHD DACC.
