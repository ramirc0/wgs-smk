# wgs-smk

A [Snakemake](https://github.com/snakemake/snakemake) pipeline for germline
whole genome sequencing of single samples. It follows WARP's DRAGEN maximum
quality mode. Paired-end FASTQs are split into read groups, aligned with
DRAGMAP and duplicate-marked. GATK HaplotypeCaller then calls variants in
DRAGEN mode. Each sample ends as a reblocked gVCF, the input of
[joint-genotyping-smk](https://github.com/christian728504/joint-genotyping-smk.git).

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
    common.smk           config, sample sheet, path constants, helpers
    reference.smk        hash table, STR table, interval lists, dbSNP, scatter
    align.smk            read group uBAMs, DRAGMAP + MergeBamAlignment, MarkDuplicates, sort
    cram.smk             CRAM, prevalidation, ValidateSamFile
    call.smk             DRAGstr, HaplotypeCaller, hard filter, gather, reblock
    qc.smk               Picard metrics, VerifyBamID, fingerprints
  scripts/               argparse scripts, each with build_parser() + main()
  envs/pixi.toml         rule env (pixi workspace)
.tests/
  unit/                  pytest tests of the scripts
  plumbing/              1% of ENCSR410ALS: config + sample sheet
  integration/           ENCODE end-to-end test: config, sample sheet, checks
profiles/
  local/config.yaml      local execution
  slurm/config.yaml      SLURM (per-rule resources, retried with more memory)
pixi.toml                launcher env, lint and test tasks
pyproject.toml           ruff, numpydoc and pytest config
```

## Running

```bash
cp config/config.yaml.template config/config.yaml  # then edit
cp config/samples.example.tsv config/samples.tsv   # then edit
pixi shell                                         # launcher env
snakemake -n -p --profile profiles/local           # dry run
snakemake --profile profiles/local                 # local
snakemake --profile profiles/slurm                 # SLURM
```

The SLURM profile holds Weng Lab cluster settings, such as partitions, excluded
nodes and the `cascadelake` constraint. Adapt them for another cluster with the
[SLURM executor docs](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html).
Tune `set-resources` from the `benchmarks/<run_id>/` TSVs after a first run.

## Tests

`pixi run test` runs the unit tests. `pixi run lint` checks the scripts with
ruff and `numpydoc lint`.

`.tests/plumbing/` runs every 100th read pair of ENCSR410ALS, 4.8 M pairs over
its 4 NextSeq lanes. Its config lowers the large read group
cutoff so that one read group is split and gathered.

```bash
snakemake --profile profiles/slurm --configfile .tests/plumbing/config.yaml
```

`.tests/integration/` runs three ENCODE whole genome samples end to end on full
GRCh38. Two are EN-TEx donors sequenced on a NextSeq 500. The third is A549 on a
HiSeq X Ten. Its inputs live on the Weng Lab cluster, so it runs only there.

```bash
snakemake --profile profiles/slurm --configfile .tests/integration/config.yaml
pixi run --manifest-path workflow/envs/pixi.toml python .tests/integration/check_run.py
pixi run --manifest-path workflow/envs/pixi.toml python .tests/integration/compare_encode.py
```

Outputs go under `results/encode-integration/`. `check_run.py` checks WARP
invariants: the CRAM holds every read in its flowcell and lane's read group,
HaplotypeCaller ran with contamination 0, the WGS metrics exist and
ValidateSamFile found no errors. `compare_encode.py` compares the EN-TEx
genotypes with ENCODE's GRCh38 variant calls for the same donors: `ENCFF501WAL`
(ENCSR410ALS) and `ENCFF568LEP` (ENCSR456SNK).

## Inputs

The sample sheet is a TSV with a header row. `sample_id` names the sample. `r1`
and `r2` are its gzipped paired-end FASTQs. Lines that start with `#` are
comments. Each `sample_id` MUST be unique. Every read name MUST hold 7
(Casava 1.8) or 8 (with a UMI) colon-separated fields. Fields 3 and 4, the
flowcell and lane, set the read group.

## References

Configs read the references from `/data/zusers/ramirezc/static/hg38/`. `B` is
`https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0`.

| Config key | Static path | Download |
|---|---|---|
| `fasta` | `genome/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta`, with `.fai` and `.dict` | [ENCODE](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz), gunzipped |
| `haplotype_database_file` | `broad-references/v0/Homo_sapiens_assembly38.haplotype_database.txt` | `B/` same name |
| `dbsnp_vcf` | `broad-references/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz` | `B/` same name |
| `calling_interval_list` | `broad-references/v0/wgs_calling_regions.hg38.interval_list` | `B/` same name |
| `evaluation_interval_list` | `broad-references/v0/wgs_evaluation_regions.hg38.interval_list` | `B/` same name |
| `coverage_interval_list` | `broad-references/v0/wgs_coverage_regions.hg38.interval_list` | `B/` same name |
| `contamination_sites` | `broad-references/v0/contamination-resources/1000g/1000g.phase3.100k.b38.vcf.gz.dat`, the prefix of `.UD`, `.bed` and `.mu` | `B/contamination-resources/1000g/` |

The pipeline builds the files that depend on the FASTA once, under
`results/refs/`:

- `dragen_reference/`: the DRAGEN hash table.
- `str_table.tsv`: the STR table, from `gatk ComposeSTRTableFile`.
- `coverage.interval_list` and `evaluation.interval_list`: the Broad lists with
  the FASTA's dictionary as header. Picard checks it against the BAM's.
- `dbsnp.vcf.gz`: dbSNP with the FASTA's contigs in its header.

## Configuration

| Key | Meaning |
|---|---|
| `samples` | Sample sheet path. |
| `outdir` | Results root. |
| `run_id` | Run name. Outputs go under `results/<run_id>/`. |
| `pixi_env` | Environment of `workflow/envs/pixi.toml` that the rules use. |
| `picard_tmp_dir` | JVM temporary directory. GKL and snappy load native libraries from it, so it MUST NOT sit on a `noexec` mount, as `/tmp` is on the compute nodes. Sort spills go to `$TMPDIR` instead. |
| `references.*` | See [References](#references). |
| `platform` | Read group platform. Default `Generic`. |
| `compression_level` | BAM compression level. Default `2`. |
| `lod_threshold` | CrosscheckFingerprints LOD threshold. Default `-20.0`. |
| `cross_check_fingerprints_by` | CrosscheckFingerprints `CROSSCHECK_BY`. Default `READGROUP`. |
| `contamination_underestimation_factor` | Divisor of the FREEMIX that `check_contamination` logs. Default `0.75`. |
| `split_large_read_group.cutoff_gib` | Split read group uBAMs larger than this. WARP uses 20. |
| `split_large_read_group.reads_per_file` | Reads per split uBAM. WARP uses 48 M. |
| `scatter_interval_list.scatter_count` | Number of HaplotypeCaller scatters. WARP uses 50. |
| `scatter_interval_list.break_bands_at_multiples_of` | IntervalListTools `BREAK_BANDS_AT_MULTIPLES_OF`. WARP uses 1000000. |

## Steps

References, in `reference.smk`:

- `reference_index` writes the FASTA index and sequence dictionary when they are missing.
- `dragen_index` builds the DRAGEN hash table.
- `compose_str_table` builds the STR table.
- `reheader_interval_list` gives the coverage and evaluation lists the FASTA's dictionary.
- `prepare_dbsnp` gives dbSNP the FASTA's contigs.
- `scatter_interval_list` splits the calling regions with IntervalListTools.

Alignment, in `align.smk`:

- `unmapped_bam` routes each read pair to its read group, `<flowcell>.<lane>`,
  and runs FastqToSam per read group. It splits a read group uBAM larger than
  `cutoff_gib` into files of `reads_per_file` reads, as WARP's
  SplitLargeReadGroup does.
- `read_groups` keeps the table of read groups and their split uBAMs. It is the
  only checkpoint.
- `unmapped_bam_to_aligned` streams DRAGMAP's alignments of one uBAM into MergeBamAlignment.
- `gather_unsorted_bam_files` joins the alignments of a split read group.
- `mark_duplicates` runs MarkDuplicates over all read groups.
- `sort_bam` coordinate-sorts and indexes it. This is the final BAM.

CRAM, in `cram.smk`:

- `convert_to_cram` writes a CRAM 3.0 file with its index and md5.
- `check_prevalidation` extracts the duplication and chimerism rates. Above
  0.30 or 0.15, the sample is outlier data.
- `validate_sam_file` runs ValidateSamFile. It skips mate validation for outlier data.

Calling, in `call.smk`:

- `calibrate_dragen_str_model` fits the DRAGstr model.
- `haplotype_caller` calls each scatter in DRAGEN mode with contamination 0.
- `dragen_hardfilter_gvcf` flags `QUAL < 10.4139` as `DRAGENHardQUAL`.
- `merge_gvcfs` gathers the scatters by block copy.
- `reblock` runs ReblockGVCF.
- `validate_gvcf` runs ValidateVariants and touches `<sample>.OK`.

QC, in `qc.smk`:

- `collect_quality_yield_metrics` runs on each read group uBAM.
- `collect_unsorted_readgroup_bam_quality_metrics` runs on each read group's alignments.
- `cross_check_fingerprints` checks that all read groups come from one individual.
- `check_contamination` estimates contamination with VerifyBamID.
- `collect_readgroup_bam_quality_metrics`, `collect_aggregation_metrics`,
  `collect_wgs_metrics`, `collect_raw_wgs_metrics` and
  `calculate_readgroup_checksum` run on the final BAM.
- `collect_variant_calling_metrics` runs on the reblocked gVCF.

## Outputs

Each sample writes to `results/<run_id>/<sample>/<step>/`. Read group steps
name their files `<flowcell>.<lane>`. The files built from the references go to
`results/refs/`. Logs and benchmarks go to `logs/<run_id>/<step>/` and
`benchmarks/<run_id>/<step>/`. To list one step across samples, glob it:
`results/<run_id>/*/reblock/*.rb.g.vcf.gz`.

The run keeps the final BAM, the CRAM, the reblocked gVCF, the DRAGstr model and
every metrics file. The steps marked temporary below are deleted once no job
needs them. WARP keeps the same set without the BAM and the DRAGstr model.

| Step directory | Key outputs |
|---|---|
| `unmapped_bam/` | uBAM per read group, and split uBAMs of large read groups (temporary) |
| `read_groups/` | Read groups with their barcode, read pairs and split uBAMs |
| `collect_quality_yield_metrics/` | Quality yield metrics per read group |
| `unmapped_bam_to_aligned/` | Aligned BAM per uBAM, merged with it (temporary) |
| `gather_unsorted_bam_files/` | Aligned BAM per split read group (temporary) |
| `collect_unsorted_readgroup_bam_quality_metrics/` | Base distribution, insert size, quality by cycle and quality distribution per read group |
| `mark_duplicates/` | Duplicate-marked BAM (temporary), duplication metrics |
| `sort_bam/` | Final BAM, index, md5 |
| `cross_check_fingerprints/` | Fingerprint metrics |
| `check_contamination/` | VerifyBamID `.selfSM` |
| `collect_readgroup_bam_quality_metrics/` | Read group alignment summary and GC bias metrics |
| `collect_aggregation_metrics/` | Alignment, insert size, artifact, GC bias and quality distribution metrics |
| `collect_wgs_metrics/`, `collect_raw_wgs_metrics/` | WGS and raw WGS metrics over the coverage regions |
| `calculate_readgroup_checksum/` | Read group checksum |
| `convert_to_cram/` | CRAM, its index and md5 |
| `check_prevalidation/` | Duplication and chimerism rates, outlier flag |
| `validate_sam_file/` | ValidateSamFile report |
| `calibrate_dragen_str_model/` | DRAGstr model |
| `haplotype_caller/` | Per-scatter gVCFs (temporary) |
| `dragen_hardfilter_gvcf/` | Hard-filtered per-scatter gVCFs (temporary) |
| `merge_gvcfs/` | Gathered gVCF (temporary) |
| `reblock/` | Reblocked gVCF (`.rb.g.vcf.gz`), ready for joint genotyping |
| `validate_gvcf/` | `<sample>.OK` once ValidateVariants passes |
| `collect_variant_calling_metrics/` | Variant calling detail and summary metrics |

## Deviations from WARP

The pipeline matches WARP `WholeGenomeGermlineSingleSample` 3.3.8 in
`dragen_maximum_quality_mode` except here:

- The reference is ENCODE's `GRCh38_no_alt_analysis_set` instead of
  `Homo_sapiens_assembly38`, which has ALT, HLA and decoy contigs.
- Input is a FASTQ pair per sample. WARP takes one uBAM per read group.
- `samtools sort` replaces Picard SortSam. Its stable sort of name-ordered
  input keeps htsjdk's order of tied records.
- ReblockGVCF runs with `--add-site-filters-to-genotype`. Without it, ReblockGVCF
  strips `DRAGENHardQUAL` from every variant. With it, the filter survives as
  the genotype's FT. GenotypeGVCFs ignores both.
- `QUAL < 10.4139` also flags every reference block, because a missing QUAL
  counts as -10. WARP does the same. GenotypeGVCFs ignores site filters.
- ValidateSamFile also ignores `MATE_NOT_FOUND`.
- MarkDuplicates parses tile, x and y from fields 5-7 of 7- and 8-field read
  names. WARP relies on Picard's default parser, which takes no 8-field name.
- The tools are newer. WARP's images run GATK 4.6.1.0, Picard 2.26.10,
  samtools 1.11 and DRAGMAP 1.2.1.

## Environment

The launcher env (`pixi.toml` at the root) has Snakemake main and its plugins
for pixi deployment and SLURM. The rule env lives in a second pixi workspace,
`workflow/envs/`. Its `pixi.lock` pins every package: GATK 4.7.0.0, Picard
3.5.0, samtools, bcftools and htslib 1.24, DRAGMAP 1.3.0 and VerifyBamID 2.0.3
from bioconda, on OpenJDK 25. pixi MUST be on `PATH`. An `onstart` hook
installs the env before the first job. Parallel SLURM jobs would otherwise race
to install it.

Bioconda's `gatk4-main` pins OpenJDK 17, which slows Picard by about 30%.
`workflow/envs/gatk4-main/` repackages the same GATK release with `openjdk >=17`.
GATK's build accepts Java 17 or later but tests only 17. On a 5 Mb region,
HaplotypeCaller wrote identical records on 17 and 25.

DRAGMAP's repository is archived. Its last commit is from 2022. WARP's DRAGEN
mode has no other mapper.

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
- `input:` and `output:` are literal path strings, or functions of the
  `read_groups` checkpoint. All paths are relative.
- Java rules take their JVM options from `java_options()` in `common.smk`.
- Rules tee logs to file and terminal (`exec &> >(tee {log:q})`).
- Every rule has `log:` and `benchmark:`. Each names its pixi env in `software:`.
- Shell blocks end with `sync {output:q}`. On a shared filesystem, Snakemake's
  first stat of an output from the submit host otherwise waits for the
  writing node to flush it.
- Docstrings follow numpydoc. Neither ruff nor `numpydoc lint` parses `.smk` files.

## Citations

Behera, Sairam, Severine Catreux, Massimiliano Rossi, et al. 2024. “Comprehensive Genome Analysis and Variant Detection at Scale Using DRAGEN.” Nature Biotechnology 43 (7): 1177–91. https://doi.org/10.1038/s41587-024-02382-1.

Danecek, Petr, James K. Bonfield, Jennifer Liddle, et al. 2021. “Twelve Years of SAMtools and BCFtools.” GigaScience 10 (2): giab008. https://doi.org/10.1093/gigascience/giab008.

McKenna, Aaron, Matthew Hanna, Eric Banks, et al. 2010. “The Genome Analysis Toolkit: A MapReduce Framework for Analyzing next-Generation DNA Sequencing Data.” Genome Research 20 (9): 1297–303. https://doi.org/10.1101/gr.107524.110.

Mölder, Felix, Kim Philipp Jablonski, Brice Letcher, et al. 2021. “Sustainable Data Analysis with Snakemake.” F1000Research 10 (April): 33. https://doi.org/10.12688/f1000research.29032.2.

Warp: Warp Analysis Research Pipelines. n.d. Github. Accessed October 10, 2026. https://github.com/broadinstitute/warp.

Zhang, Fan, Matthew Flickinger, Sarah A. Gagliano Taliun, et al. 2020. “Ancestry-Agnostic Estimation of DNA Sample Contamination from Sequence Reads.” Genome Research 30 (2): 185–94. https://doi.org/10.1101/gr.246934.118.

## Questions

If you have any questions or would like to provide constructive feedback, please open an issue or reach out to the MOHD DACC.
