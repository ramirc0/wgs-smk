# CLAUDE.md

Germline WGS Snakemake workflow on the snakemake-template layout.
README.md holds the layout, steps, outputs, and style conventions. Follow it.

## Commands

```bash
pixi run lint                                # ruff check, ruff format --check, numpydoc lint
pixi install --all --manifest-path workflow/envs/pixi.toml  # once per checkout; compiles the source packages
pixi run snakemake -n -p --profile profiles/local
pixi run snakemake --profile profiles/slurm
```

## Two pixi workspaces

- Root `pixi.toml` is the launcher: Snakemake and its plugins, plus the `lint` env.
- `workflow/envs/pixi.toml` holds the rule envs `default` and `gatk`. Rule dependencies go there, never in the root.
- `pyproject.toml` is ruff and numpydoc config only. pixi ignores it.
- Rules use `locked=True`. After a manifest or recipe edit, run `pixi lock` in `workflow/envs/`.
- Any manifest or lock change reruns every rule.
- `SOFTWARE_ENV` and `GATK_ENV` resolve `../envs` relative to the rule file. New rule files MUST sit in `workflow/rules/`.

## Source packages

- `workflow/envs/dragmap/` and `workflow/envs/verifybamid2/` are pixi-build recipes (rattler-build backend).
- They build DRAGMAP `c4cb6c94` and VerifyBamID `4b362f51` with GCC 9, the image's CFLAGS/CXXFLAGS/LDFLAGS and Ubuntu GCC's default hardening flags.
- Conda's default `CPPFLAGS` MUST stay unset. `-DNDEBUG` compiles out the `posix_memalign` that DRAGMAP wraps in `BOOST_ASSERT`.
- Without a prior `pixi install --all`, the first job compiles both.

## Outputs must not drift

- Pins match `clarity001/wgs-smk:latest` pulled 2025-12-22 and `us.gcr.io/broad-gatk/gatk:latest` pulled 2025-07-31. Their SIFs are in `/data/zusers/ramirezc/wgs-smk/workflow/rules/envs/`. Docker Hub `:latest` and the old repo Dockerfile do not match them.
- `picard.jar` is byte-identical to the image's. bioconda's GATK jar has the image's classes compiled by a newer javac.
- Tag `pre-template` is the last commit before the template port.
- README.md "Known issues" lists legacy bugs kept on purpose. A fix changes outputs and needs the user's approval.

## Gotchas

- Inside a double-quoted shell string, `:q` embeds literal quotes. GATK `--java-options` and MergeBamAlignment's `PROGRAM_GROUP_COMMAND_LINE` interpolate without it.
- `-Xlog:gc*:file={log}` makes the JVM rotate the tee'd log. The rule output ends in `{log}.0`, the GC log in `{log}`.
- Picard scratch goes to `picard_tmp_dir`, never `/tmp`. Picard fails on a `noexec` `/tmp`.
- `* attempt` in `profiles/slurm/config.yaml` is intentional. JVM heaps follow `mem_mb`, so a retry also doubles them.
- `set-threads` sets `OMP_NUM_THREADS` to the job's threads. Before the port it was 1.
- `reference_index` takes the FASTA as `ancient()`. Existing index files stay even when the FASTA is newer. A regenerated `.dict` would change BAM `@SQ` lines.
- Resources MUST be flat values. Never derive `mem_mb` or `runtime` from input size.
- A dry run with the SLURM profile defaults to `--cores 1` and shows every `{threads}` as 1. Pass `--cores 104`.
- CLI `--set-threads` replaces the profile's whole mapping. Pass every rule.

## Verifying a change

No test suite. Diff against `pre-template` on a real sample:

1. Take the first 10M read pairs of EG100001.
2. Run `git archive pre-template` with the two production SIFs, linked in a scratch `apptainer-prefix` under the md5 of their image URLs.
3. Run this workflow with the same per-rule threads.
4. Compare `samtools view` records and headers without `@PG`, `bcftools view` records and headers without `##GATKCommandLine`, and Picard metrics without `#` lines. BQSR tables, `.selfSM`, interval lists, `.dragstr` and the validation report MUST match as text.
5. md5 files and PDFs differ: they embed paths and dates.
