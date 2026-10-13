"""Extract the duplication and chimerism rates that set ValidateSamFile's mode."""

import argparse
import csv
import io
import subprocess
from pathlib import Path

# WARP's BamToCram: ValidateSamFile skips mate validation above these, where
# it can run out of memory.
MAX_DUPLICATION = 0.30
MAX_CHIMERISM = 0.15


def build_parser():
    """Build the check_prevalidation.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the check_prevalidation.py flags.
    """
    p = argparse.ArgumentParser(description="Extract duplication and chimerism.")
    p.add_argument(
        "--duplication-metrics", required=True, help="MarkDuplicates metrics."
    )
    p.add_argument(
        "--alignment-summary-metrics",
        required=True,
        help="CollectMultipleMetrics alignment summary metrics.",
    )
    p.add_argument(
        "--duplication-value", required=True, help="Text file for PERCENT_DUPLICATION."
    )
    p.add_argument(
        "--chimerism-value", required=True, help="Text file for PCT_CHIMERAS."
    )
    p.add_argument(
        "--is-outlier-data",
        required=True,
        help="Text file for `true` or `false`: either rate is above its maximum.",
    )
    return p


def write_last_value(table, column, path):
    """Write `column` of the last row of a tab-separated table, without a newline.

    Parameters
    ----------
    table : str
        Tab-separated text with a header row.
    column : str
        Column to write.
    path : str
        File to write. It is not created when the table has no data rows.
    """
    for row in csv.DictReader(io.StringIO(table), delimiter="\t"):
        with open(path, "w") as f:
            f.write(row[column])


def main(argv=None):
    """Write PERCENT_DUPLICATION, the PAIR row's PCT_CHIMERAS and the outlier flag.

    `grep` selects the rows and fails when a metric header is missing.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)

    duplication = subprocess.run(
        ["grep", "-A", "1", "PERCENT_DUPLICATION", args.duplication_metrics],
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    ).stdout
    # Header plus the FIRST_OF_PAIR, SECOND_OF_PAIR and PAIR rows. Like the
    # `grep -A 3 | grep -v OF_PAIR` pipeline, only the second grep's status counts.
    chimerism_rows = subprocess.run(
        ["grep", "-A", "3", "PCT_CHIMERAS", args.alignment_summary_metrics],
        check=False,
        stdout=subprocess.PIPE,
        text=True,
    ).stdout
    chimerism = subprocess.run(
        ["grep", "-v", "OF_PAIR"],
        input=chimerism_rows,
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    ).stdout

    write_last_value(duplication, "PERCENT_DUPLICATION", args.duplication_value)
    write_last_value(chimerism, "PCT_CHIMERAS", args.chimerism_value)
    is_outlier = (
        float(Path(args.duplication_value).read_text()) > MAX_DUPLICATION
        or float(Path(args.chimerism_value).read_text()) > MAX_CHIMERISM
    )
    Path(args.is_outlier_data).write_text(str(is_outlier).lower())


if __name__ == "__main__":
    main()
