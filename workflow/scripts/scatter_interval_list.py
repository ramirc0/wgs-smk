"""Scatter an interval list with Picard IntervalListTools and number the pieces."""

import argparse
import glob
import os
import shlex
import subprocess
import sys


def build_parser():
    """Build the scatter_interval_list.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the scatter_interval_list.py flags.
    """
    p = argparse.ArgumentParser(description="Scatter an interval list.")
    p.add_argument("--java-options", required=True, help="JVM options, one string.")
    p.add_argument("--scatter-count", type=int, required=True, help="Number of pieces.")
    p.add_argument(
        "--break-bands-at-multiples-of",
        type=int,
        required=True,
        help="Break bands at multiples of this many bases.",
    )
    p.add_argument("--interval-list", required=True, help="Interval list to scatter.")
    p.add_argument("--outdir", required=True, help="Output directory.")
    return p


def main(argv=None):
    """Scatter the interval list, then prefix each piece with its 1-based index.

    IntervalListTools writes `temp_<n>_of_<N>/scattered.interval_list`. Each
    file becomes `temp_<n>_of_<N>/<n>scattered.interval_list`.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.

    Returns
    -------
    int
        Exit code: 0 on success, 1 if IntervalListTools failed.
    """
    args = build_parser().parse_args(argv)

    cmd = [
        "picard",
        *shlex.split(args.java_options),
        "IntervalListTools",
        "-SCATTER_COUNT",
        str(args.scatter_count),
        "-SUBDIVISION_MODE",
        "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW",
        "-UNIQUE",
        "true",
        "-SORT",
        "true",
        "-BREAK_BANDS_AT_MULTIPLES_OF",
        str(args.break_bands_at_multiples_of),
        "-INPUT",
        args.interval_list,
        "-OUTPUT",
        args.outdir,
    ]
    result = 0
    try:
        subprocess.check_call(cmd)
    except (OSError, subprocess.CalledProcessError) as e:
        result = 1
        print(e)

    intervals = sorted(glob.glob(f"{args.outdir}/*/*.interval_list"))
    for i, interval in enumerate(intervals):
        directory, filename = os.path.split(interval)
        os.rename(interval, os.path.join(directory, str(i + 1) + filename))
    print(len(intervals))

    if result != 0:
        print(f"An error occurred. Exiting with code {result}")
    return result


if __name__ == "__main__":
    sys.exit(main())
