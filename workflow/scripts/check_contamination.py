"""Estimate cross-sample contamination with VerifyBamID and check its .selfSM."""

import argparse
import csv
import shlex
import subprocess
import sys
from pathlib import Path


def build_parser():
    """Build the check_contamination.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the check_contamination.py flags.
    """
    p = argparse.ArgumentParser(description="Run VerifyBamID on a sorted BAM.")
    p.add_argument("--threads", type=int, required=True, help="VerifyBamID threads.")
    p.add_argument("--output", required=True, help="The .selfSM file to write.")
    p.add_argument("--bam", required=True, help="Coordinate-sorted, indexed BAM.")
    p.add_argument("--reference", required=True, help="Reference FASTA.")
    p.add_argument(
        "--svd-prefix",
        required=True,
        help="Prefix of the VerifyBamID .UD, .bed and .mu resources.",
    )
    p.add_argument(
        "--underestimation-factor",
        type=float,
        required=True,
        help="FREEMIX is divided by this factor before it is logged.",
    )
    return p


def check_self_sm(path, underestimation_factor):
    """Check the VerifyBamID likelihoods and log the scaled FREEMIX.

    Parameters
    ----------
    path : str
        The .selfSM file.
    underestimation_factor : float
        FREEMIX is divided by this factor before it is logged.

    Returns
    -------
    int
        0 for one row with nonzero likelihoods, else 1.
    """
    with Path(path).open("r") as self_sm:
        rows = list(csv.DictReader(self_sm, delimiter="\t"))
    if len(rows) != 1:
        print(f"Found {len(rows)} rows in {path}. Expected exactly 1.", file=sys.stderr)
        return 1
    row = rows[0]
    if float(row["FREELK0"]) == 0 and float(row["FREELK1"]) == 0:
        # Zero likelihoods mean no data: a very shallow BAM, or a BAM aligned
        # to another reference build than the VCF.
        print(
            "Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).",
            file=sys.stderr,
        )
        return 1
    print(float(row["FREEMIX"]) / underestimation_factor, file=sys.stderr)
    return 0


def main(argv=None):
    """Run VerifyBamID, then check its output.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.

    Returns
    -------
    int
        Exit code: 0 on success, 1 if the .selfSM check fails.
    """
    args = build_parser().parse_args(argv)

    cmd = [
        # bioconda's wrapper passes arguments through only with --SVDPrefix.
        "verifybamid2",
        "--Verbose",
        "--NumPC",
        "4",
        "--NumThread",
        str(args.threads),
        "--Output",
        args.output.removesuffix(".selfSM"),
        "--BamFile",
        args.bam,
        "--Reference",
        args.reference,
        "--SVDPrefix",
        args.svd_prefix,
    ]
    print(f"Running command:\n\n{shlex.join(cmd)}", flush=True)
    subprocess.check_call(cmd)

    print("Verifying output")
    result = check_self_sm(args.output, args.underestimation_factor)
    if result != 0:
        print(f"An error occurred. Exiting with code {result}")
    return result


if __name__ == "__main__":
    sys.exit(main())
