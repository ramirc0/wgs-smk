# wgs-smk

A [Snakemake](https://github.com/snakemake/snakemake) pipeline for germline
whole genome sequencing of single samples. Paired-end FASTQs are aligned with
DRAGMAP, duplicate-marked and base-quality recalibrated. GATK HaplotypeCaller
then calls variants in DRAGEN mode. Each sample ends as a reblocked gVCF, the
input of [joint-genotyping-smk](https://github.com/christian728504/joint-genotyping-smk.git).

> [!NOTE]
> This pipeline is based on the [WARP Whole Genome Germline Single Sample pipeline](https://broadinstitute.github.io/warp/docs/Pipelines/Whole_Genome_Germline_Single_Sample/README) by the Broad Institute, [originally licensed](https://github.com/broadinstitute/warp/blob/develop/LICENSE) under BSD. See [THIRD_PARTY_LICENSES](THIRD_PARTY_LICENSES) for details.

## Layout

```
config/
  config.yaml.template   parameters + paths; copy to config.yaml (gitignored)
  samples.example.tsv    one row per sample: sample_id, r1, r2
workflow/
  Snakefile              includes the rule files, defines rule all
  rules/
    common.smk           config, sample sheet, path constants, repeat_flag()
    reference.smk        FASTA indexes, DRAGEN hash table, interval scatter
    align.smk            uBAM, DRAGMAP + MergeBamAlignment, MarkDuplicates, sort
    bqsr.smk             BaseRecalibrator, ApplyBQSR, gathers
    cram.smk             CRAM, prevalidation, ValidateSamFile
    call.smk             DRAGstr, HaplotypeCaller, hard filter, merge, reblock
    qc.smk               Picard metrics, VerifyBamID, fingerprints
  scripts/               argparse scripts, each with build_parser() + main()
  envs/
    pixi.toml            rule envs (pixi workspace): default and gatk
    dragmap/             source package: DRAGMAP at the production commit
    verifybamid2/        source package: VerifyBamID at the production commit
.tests/integration/      ENCODE end-to-end test: config + sample sheet
profiles/
  local/config.yaml      local execution
  slurm/config.yaml      SLURM (per-rule resources, retried with more memory)
pixi.toml                launcher env + lint task
pyproject.toml           ruff + numpydoc config
```

## Running

```bash
cp config/config.yaml.template config/config.yaml  # then edit
cp config/samples.example.tsv config/samples.tsv   # then edit
pixi install --all --manifest-path workflow/envs/pixi.toml  # once; compiles DRAGMAP and VerifyBamID
pixi shell                                         # launcher env
snakemake -n -p --profile profiles/local           # dry run
snakemake --profile profiles/local                 # local
snakemake --profile profiles/slurm                 # SLURM
```

The SLURM profile holds Weng Lab cluster settings, such as partitions, excluded
nodes and the `cascadelake` constraint. Adapt them for another cluster with the
[SLURM executor docs](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html).
Tune `set-resources` from the `benchmarks/<run_id>/` TSVs after a first run.

## Integration test

`.tests/integration/` runs three ENCODE whole genome samples end to end on full
GRCh38. Two are EN-TEx donors sequenced on a NextSeq 500. The third is A549 on a
HiSeq X Ten. Its inputs live on the Weng Lab cluster, so it runs only there.

```bash
snakemake --profile profiles/slurm --configfile .tests/integration/config.yaml
```

Outputs go under `results/encode-integration/`. Compare the EN-TEx genotypes with
ENCODE's GRCh38 variant calls for the same donors: `ENCFF501WAL` (ENCSR410ALS)
and `ENCFF568LEP` (ENCSR456SNK).

## Inputs

The sample sheet is a TSV with a header row. `sample_id` names the sample. `r1`
and `r2` are its gzipped paired-end FASTQs. Lines that start with `#` are
comments. Each `sample_id` MUST be unique. The first R1 read name MUST hold 7
(Casava 1.8) or 8 (with a UMI) colon-separated fields. They set the read group
of every read.

Download the references into `resources/` (gitignored):

```bash
curl -fsSL https://users.wenglab.org/ramirezc/mohd/resources/wgs-smk.tar.zst | tar --use-compress-program=zstd -xf -
```

`reference_index` writes the `.fai`, `.dict` and `.chromsizes` files beside the
FASTA when they are missing.

## Configuration

| Key | Meaning |
|---|---|
| `samples` | Sample sheet path. |
| `outdir` | Results root. |
| `run_id` | Run name. Outputs go under `results/<run_id>/`. |
| `pixi_env` | Environment of `workflow/envs/pixi.toml` that non-GATK rules use. |
| `tmp_dir` | Scratch for the aligned BAM, the CRAM reference cache and the dbSNP copy. |
| `picard_tmp_dir` | Picard scratch. It MUST NOT sit on a `noexec` mount, as `/tmp` is on some hosts. |
| `references.fasta` | Reference genome FASTA. |
| `references.haplotype_database_file` | Haplotype map for CrosscheckFingerprints. |
| `references.known_indel_sites_vcfs` | Known indel VCFs for BaseRecalibrator. |
| `references.dbsnp_vcf` | dbSNP VCF for BaseRecalibrator, ValidateVariants and variant calling metrics. |
| `references.calling_interval_list` | Calling regions to scatter. |
| `references.evaluation_interval_list` | Target intervals for variant calling metrics. |
| `references.sequence_grouping` | One BQSR group per line, tab-separated contigs. |
| `references.str_table_file` | STR table for CalibrateDragstrModel. |
| `references.contamination_sites` | VerifyBamID `.UD`, `.bed` and `.mu` resources. |
| `platform` | Read group platform. Default `Generic`. |
| `compression_level` | BAM compression level. Default `2`. |
| `unmap_contaminant_reads` | MergeBamAlignment `UNMAP_CONTAMINANT_READS`. Default `false`. |
| `lod_threshold` | CrosscheckFingerprints LOD threshold. Default `-20.0`. |
| `cross_check_fingerprints_by` | CrosscheckFingerprints `CROSSCHECK_BY`. Default `READGROUP`. |
| `contamination_underestimation_factor` | Divisor of the FREEMIX that `check_contamination` logs. Default `0.75`. |
| `bin_base_qualities` | Quantize recalibrated qualities to 10, 20 and 30. |
| `somatic` | With `bin_base_qualities`, also quantize to 30, 40 and 50. |
| `scatter_interval_list.scatter_count` | Number of HaplotypeCaller scatters. |
| `scatter_interval_list.break_bands_at_multiples_of` | IntervalListTools `BREAK_BANDS_AT_MULTIPLES_OF`. |

## Steps

References, in `reference.smk`:

- `reference_index` writes the FASTA index, sequence dictionary and chromosome sizes.
- `dragen_index` builds the DRAGEN hash table.
- `scatter_interval_list` splits the calling regions with IntervalListTools.

Alignment, in `align.smk`:

- `unmapped_bam` converts the FASTQs to a uBAM with FastqToSam.
- `unmapped_bam_to_aligned` aligns the uBAM with DRAGMAP and merges the alignments back with MergeBamAlignment.
- `mark_duplicates` runs MarkDuplicates on the queryname-grouped BAM.
- `sort_bam` coordinate-sorts and indexes it.

BQSR, in `bqsr.smk`:

- `base_recalibration` runs BaseRecalibrator once per sequence group.
- `gather_bqsr_reports` merges the group reports.
- `apply_bqsr` recalibrates each sequence group.
- `gather_recalibrated_bams` concatenates the groups into the final BAM.

CRAM, in `cram.smk`:

- `convert_to_cram` writes a CRAM 3.0 file, its md5 and its index.
- `check_prevalidation` extracts the duplication and chimerism rates.
- `validate_sam_file` runs ValidateSamFile. It skips mate validation when
  duplication exceeds 0.30 or chimerism exceeds 0.15.

Calling, in `call.smk`:

- `calibrate_dragen_str_model` fits the DRAGstr model.
- `haplotype_caller` calls each scatter in DRAGEN mode. VerifyBamID's FREEMIX is its `-contamination`.
- `dragen_hardfilter_gvcf` flags `QUAL < 10.4139` as `DRAGENHardQUAL`.
- `merge_gvcfs` merges the scatters.
- `sort_bamout` and `merge_bamout` sort and merge the HaplotypeCaller bamouts.
- `reblock` runs ReblockGVCF.
- `validate_gvcf` runs ValidateVariants and touches `<sample>.OK`.

QC, in `qc.smk`:

- `collect_quality_yield_metrics` runs on the uBAM.
- `collect_unsorted_readgroup_bam_quality_metrics` runs on the merged, unsorted BAM.
- `cross_check_fingerprints` checks that all read groups come from one individual.
- `check_contamination` estimates contamination with VerifyBamID.
- `collect_readgroup_bam_quality_metrics`, `collect_aggregation_metrics` and
  `calculate_readgroup_checksum` run on the final BAM.
- `collect_variant_calling_metrics` runs on the reblocked gVCF.

## Outputs

Each sample writes to `results/<run_id>/<sample>/<step>/`. The interval scatter
depends on the run config, so it goes to `results/<run_id>/scatter_interval_list/`.
The DRAGEN hash table does not depend on the run. It goes to
`results/refs/dragen_reference/`. Logs and benchmarks go to
`logs/<run_id>/<step>/` and `benchmarks/<run_id>/<step>/`, one file per sample.
To list one step across samples, glob it: `results/<run_id>/*/reblock/*.rb.g.vcf.gz`.

| Step directory | Key outputs |
|---|---|
| `unmapped_bam/` | uBAM with read group metadata |
| `collect_quality_yield_metrics/` | Quality yield metrics of the uBAM |
| `unmapped_bam_to_aligned/` | Aligned BAM merged with the uBAM |
| `collect_unsorted_readgroup_bam_quality_metrics/` | Insert size metrics and histogram PDF |
| `mark_duplicates/` | Duplicate-marked BAM, duplication metrics |
| `sort_bam/` | Coordinate-sorted BAM, index, md5 |
| `cross_check_fingerprints/` | Fingerprint metrics |
| `check_contamination/` | VerifyBamID `.selfSM` |
| `base_recalibrator/` | Per-group recalibration tables |
| `gather_bqsr_reports/` | Merged recalibration table |
| `apply_bqsr/` | Per-group recalibrated BAMs |
| `gather_recalibrated_bams/` | Final BAM |
| `collect_readgroup_bam_quality_metrics/` | Read group alignment summary and GC bias metrics |
| `collect_aggregation_metrics/` | Alignment, insert size, artifact, GC bias and quality distribution metrics |
| `calculate_readgroup_checksum/` | Read group checksum |
| `convert_to_cram/` | CRAM, its index and md5 |
| `check_prevalidation/` | Duplication and chimerism values |
| `validate_sam_file/` | ValidateSamFile report |
| `calibrate_dragen_str_model/` | DRAGstr model |
| `haplotype_caller/` | Per-scatter gVCFs and bamouts |
| `dragen_hardfilter_gvcf/` | Hard-filtered per-scatter gVCFs |
| `merge_gvcfs/` | Merged gVCF |
| `sort_bamout/`, `merge_bamout/` | Sorted per-scatter and merged bamouts |
| `reblock/` | Reblocked gVCF (`.rb.g.vcf.gz`), ready for joint genotyping |
| `validate_gvcf/` | `<sample>.OK` once ValidateVariants passes |
| `collect_variant_calling_metrics/` | Variant calling detail and summary metrics |

## Environment

The launcher env (`pixi.toml` at the root) has Snakemake main and its plugins
for pixi deployment and SLURM. The rule envs live in a second pixi workspace,
`workflow/envs/`. Its `pixi.lock` pins every package. `default` runs most rules.
`gatk` runs GATK 4.6.2.0 on OpenJDK 17. pixi MUST be on `PATH`.

The tool versions match the images that produced the production results:
`clarity001/wgs-smk:latest` as pulled on 2025-12-22 and
`us.gcr.io/broad-gatk/gatk:latest` as pulled on 2025-07-31. bioconda's DRAGMAP
and VerifyBamID differ from that image's builds in code. So `default` builds
both from source at the image's commits, with the image's compiler flags.
`pixi install --all` in `workflow/envs/` compiles them once. Without it, the
first job compiles them.

Add a dependency with `pixi add <pkg>` in `workflow/envs/`. After a hand edit,
run `pixi lock` there (rules use `locked=True`). Any change to the manifest or
lock reruns every rule. So does moving the workdir, because the env hash includes
the workspace's absolute path.

## Conventions

The workflow follows the
[Nextstrain Snakemake style guide](https://docs.nextstrain.org/en/latest/reference/snakemake-style-guide.html)
and [Snakemake's best practices](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html):

- No `run:` blocks. Custom logic lives in `workflow/scripts/`, called from `shell:`.
- Raw, triple-quoted `shell` blocks (`r"""`), one option per line.
- Every interpolation uses `:q` quoting, except inside a double-quoted string
  such as GATK's `--java-options`.
- Config reaches `shell` through `params:` (a `lambda`), never directly.
- `input:` and `output:` are literal path strings. All paths are relative.
- Rules tee logs to file and terminal (`exec &> >(tee {log:q})`).
- Every rule has `log:` and `benchmark:`. Each names its pixi env in `software:`.
- Docstrings follow numpydoc. `pixi run lint` checks the scripts with ruff and
  `numpydoc lint`. Neither tool parses `.smk` files.

## Known issues

These legacy bugs stay to keep outputs identical to the production results.
A fix changes outputs.

- `haplotype_caller` passes the raw FREEMIX from the `.selfSM` as
  `-contamination`. WARP passes FREEMIX divided by
  `contamination_underestimation_factor` (0.75). Here the factor only scales the
  value that `check_contamination.py` logs.
- `check_contamination.py` never checks that the `.selfSM` has one row. WARP's
  check sits inside its row loop. This pipeline's copy had it after the loop,
  where its counter was always 1. An empty `.selfSM` fails with a Python
  traceback.
- `check_contamination.py` takes the VerifyBamID prefix with
  `rstrip(".selfSM")`, which strips trailing characters, not the suffix. With
  the `.preBqsr.selfSM` name it stops at the `r` of `preBqsr`.
- `unmapped_bam.py` builds one read group from the first R1 header. Every read
  gets it, even when a FASTQ mixes lanes or flowcells. The library is always
  `MOHD-<sample>`.
- `unmapped_bam_to_aligned` writes the aligned BAM to
  `<tmp_dir>/<sample>.aligned.bam`. Two runs of one sample with the same
  `tmp_dir` share that path. A failed job leaves the file behind.
- `base_recalibration` declares the haplotype database as an input but never
  reads it. A change to that file reruns BQSR and everything downstream.
- `calculate_readgroup_checksum` passes `INPUT= <bam>` and `OUTPUT= <md5>`, a
  legacy Picard syntax, kept verbatim.

## References

Behera, Sairam, Severine Catreux, Massimiliano Rossi, et al. 2024. “Comprehensive Genome Analysis and Variant Detection at Scale Using DRAGEN.” Nature Biotechnology 43 (7): 1177–91. https://doi.org/10.1038/s41587-024-02382-1.

Danecek, Petr, James K. Bonfield, Jennifer Liddle, et al. 2021. “Twelve Years of SAMtools and BCFtools.” GigaScience 10 (2): giab008. https://doi.org/10.1093/gigascience/giab008.

McKenna, Aaron, Matthew Hanna, Eric Banks, et al. 2010. “The Genome Analysis Toolkit: A MapReduce Framework for Analyzing next-Generation DNA Sequencing Data.” Genome Research 20 (9): 1297–303. https://doi.org/10.1101/gr.107524.110.

Mölder, Felix, Kim Philipp Jablonski, Brice Letcher, et al. 2021. “Sustainable Data Analysis with Snakemake.” F1000Research 10 (April): 33. https://doi.org/10.12688/f1000research.29032.2.

Warp: Warp Analysis Research Pipelines. n.d. Github. Accessed October 10, 2026. https://github.com/broadinstitute/warp.

Zhang, Fan, Matthew Flickinger, Sarah A. Gagliano Taliun, et al. 2020. “Ancestry-Agnostic Estimation of DNA Sample Contamination from Sequence Reads.” Genome Research 30 (2): 185–94. https://doi.org/10.1101/gr.246934.118.

## Questions

If you have any questions or would like to provide constructive feedback, please open an issue or reach out to the MOHD DACC.
