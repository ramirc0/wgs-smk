"""Check a finished run against WARP invariants.

For each sample, the final CRAM MUST hold every read FastqToSam wrote, each in
the read group of its flowcell and lane. HaplotypeCaller MUST run with
contamination 0, as WARP does in DRAGEN mode. The WGS metrics MUST exist and
ValidateSamFile MUST find no errors. Run from the repository root:

    pixi run --manifest-path workflow/envs/pixi.toml python .tests/integration/check_run.py
"""

import argparse
import csv
import gzip
import re
import subprocess
import sys
from pathlib import Path

FASTA = (
    "/data/zusers/ramirezc/static/hg38/genome/"
    "GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
)
# Primary records only. Decode just QNAME, FLAG and aux tags (CRAM
# required_fields 0x803) and keep only RG, then count reads and RG mismatches
# in one pass.
READ_GROUP_CHECK = r"""
{
    n++
    split($1, name, ":")
    if ($12 != "RG:Z:" name[3] "." name[4]) { bad++; if (!example) example = $1 " " $12 }
}
END { print n + 0 "\t" bad + 0 "\t" example }
"""


def build_parser():
    """Build the check_run.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the check_run.py flags.
    """
    p = argparse.ArgumentParser(description="Check a run against WARP invariants.")
    p.add_argument("--run", default="results/encode-integration", help="Run directory.")
    p.add_argument("--fasta", default=FASTA, help="Reference FASTA of the CRAMs.")
    p.add_argument("--threads", type=int, default=8, help="samtools threads.")
    p.add_argument(
        "--sample",
        action="append",
        help="Sample to check. Repeat for several. Default: every sample.",
    )
    return p


def total_reads(sample_dir):
    """Sum TOTAL_READS over a sample's quality yield metrics.

    Parameters
    ----------
    sample_dir : pathlib.Path
        The sample's directory in the run.

    Returns
    -------
    int
        Reads in the unmapped BAMs, over all read groups.
    """
    total = 0
    for path in sorted(sample_dir.glob("collect_quality_yield_metrics/*.metrics")):
        lines = [
            line for line in path.read_text().splitlines() if line and line[0] != "#"
        ]
        rows = csv.DictReader(lines, delimiter="\t")
        total += int(next(rows)["TOTAL_READS"])
    return total


def scan_cram(cram, fasta, threads):
    """Count a CRAM's primary reads and those whose RG differs from their name.

    Parameters
    ----------
    cram : pathlib.Path
        Final CRAM.
    fasta : str
        Its reference FASTA.
    threads : int
        Threads for samtools.

    Returns
    -------
    tuple of (int, int, str)
        Primary reads, reads whose RG is not `<flowcell>.<lane>` of their name,
        and the first such read's name and RG.
    """
    pipeline = (
        f"samtools view -@ {threads} -F 0x900 -T {fasta}"
        f" --input-fmt-option required_fields=0x803 --keep-tag RG {cram}"
        f" | mawk -F '\\t' '{READ_GROUP_CHECK}'"
    )
    out = subprocess.run(
        ["bash", "-o", "pipefail", "-c", pipeline],
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    ).stdout
    reads, mismatched, example = out.rstrip("\n").split("\t")
    return int(reads), int(mismatched), example


def contamination(gvcf):
    """Return the contamination HaplotypeCaller ran with.

    Parameters
    ----------
    gvcf : pathlib.Path
        Reblocked gVCF. Its header keeps HaplotypeCaller's command line.

    Returns
    -------
    float or None
        `--contamination-fraction-to-filter`, or None without a HaplotypeCaller line.
    """
    with gzip.open(gvcf, "rt") as f:
        for line in f:
            if not line.startswith("##"):
                return None
            if line.startswith("##GATKCommandLine=<ID=HaplotypeCaller,"):
                match = re.search(r"contamination-fraction-to-filter (\S+)", line)
                return float(match.group(1))
    return None


def check_sample(sample_dir, fasta, threads):
    """Run every check on one sample.

    Parameters
    ----------
    sample_dir : pathlib.Path
        The sample's directory in the run.
    fasta : str
        Reference FASTA of the CRAM.
    threads : int
        Threads for samtools.

    Returns
    -------
    list of tuple of (str, bool, str)
        Check name, whether it passed, and a detail.
    """
    sample = sample_dir.name
    expected = total_reads(sample_dir)
    reads, mismatched, example = scan_cram(
        sample_dir / f"convert_to_cram/{sample}.cram", fasta, threads
    )
    hc_contamination = contamination(sample_dir / f"reblock/{sample}.rb.g.vcf.gz")
    wgs_metrics = [
        sample_dir / f"collect_wgs_metrics/{sample}.wgs_metrics",
        sample_dir / f"collect_raw_wgs_metrics/{sample}.raw_wgs_metrics",
    ]
    missing = [p.name for p in wgs_metrics if not p.is_file()]
    report = sample_dir / f"validate_sam_file/{sample}.validation_report"
    first_line = report.read_text().split("\n", 1)[0]
    return [
        ("cram_reads", reads == expected, f"{reads} in CRAM, {expected} in uBAMs"),
        ("read_groups", mismatched == 0, f"{mismatched} mismatched {example}".rstrip()),
        ("hc_contamination", hc_contamination == 0, str(hc_contamination)),
        (
            "wgs_metrics",
            not missing,
            "missing " + ", ".join(missing) if missing else "",
        ),
        ("validation", first_line == "No errors found", first_line),
    ]


def main(argv=None):
    """Print one row per sample and check.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.

    Returns
    -------
    int
        Exit code: 0 if every check passed, else 1.
    """
    args = build_parser().parse_args(argv)
    run = Path(args.run)
    samples = args.sample or sorted(
        p.parent.parent.name for p in run.glob("*/convert_to_cram/*.cram")
    )
    failed = False
    print("sample\tcheck\tresult\tdetail", flush=True)
    for sample in samples:
        for name, passed, detail in check_sample(
            run / sample, args.fasta, args.threads
        ):
            failed |= not passed
            print(
                f"{sample}\t{name}\t{'PASS' if passed else 'FAIL'}\t{detail}",
                flush=True,
            )
    return int(failed)


if __name__ == "__main__":
    sys.exit(main())
