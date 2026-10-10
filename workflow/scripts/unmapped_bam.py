"""Convert paired FASTQs to an unmapped BAM with Picard FastqToSam."""

import argparse
import gzip
import shlex
import subprocess
import sys


def build_parser():
    """Build the unmapped_bam.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the unmapped_bam.py flags.
    """
    p = argparse.ArgumentParser(description="Convert paired FASTQs to a uBAM.")
    p.add_argument("--mem-mb", type=int, required=True, help="Java heap in MB.")
    p.add_argument("--tmp-dir", required=True, help="Picard temporary directory.")
    p.add_argument("--sample", required=True, help="Sample ID.")
    p.add_argument("--platform", required=True, help="Read group platform.")
    p.add_argument("--r1", required=True, help="Gzipped R1 FASTQ.")
    p.add_argument("--r2", required=True, help="Gzipped R2 FASTQ.")
    p.add_argument("--output", required=True, help="Unmapped BAM to write.")
    return p


def fastq_to_sam_command(mem_mb, tmp_dir, sample, platform, r1, r2, output):
    """Build the FastqToSam command, with the read group from the first R1 read.

    Parameters
    ----------
    mem_mb : int
        Java heap in MB.
    tmp_dir : str
        Picard temporary directory.
    sample : str
        Sample ID.
    platform : str
        Read group platform.
    r1 : str
        Gzipped R1 FASTQ. Its first header MUST hold 8 colon-separated fields.
    r2 : str
        Gzipped R2 FASTQ.
    output : str
        Unmapped BAM to write.

    Returns
    -------
    list of str
        FastqToSam command.
    """
    with gzip.open(r1, "rt") as f:
        first_line = f.readline()
    whitespace = first_line.index(" ")
    [instrument_id, _, flowcell_id, lane_no, _, _, _, _] = (
        first_line[:whitespace].lstrip("@").split(":")
    )
    return [
        "picard",
        f"-Xms{mem_mb}m",
        f"-Xmx{mem_mb}m",
        f"-Dorg.xerial.snappy.tempdir={tmp_dir}",
        f"-Djava.io.tmpdir={tmp_dir}",
        "FastqToSam",
        "-FASTQ",
        r1,
        "-FASTQ2",
        r2,
        "-OUTPUT",
        output,
        "-READ_GROUP_NAME",
        f"{flowcell_id}.{lane_no}",
        "-SAMPLE_NAME",
        sample,
        "-LIBRARY_NAME",
        f"MOHD-{sample}",
        "-PLATFORM_UNIT",
        instrument_id,
        "-PLATFORM",
        platform,
        "-SORT_ORDER",
        "queryname",
    ]


def main(argv=None):
    """Run FastqToSam on the read pair.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.

    Returns
    -------
    int
        Exit code: 0 on success, 1 on any error.
    """
    args = build_parser().parse_args(argv)
    try:
        cmd = fastq_to_sam_command(
            args.mem_mb,
            args.tmp_dir,
            args.sample,
            args.platform,
            args.r1,
            args.r2,
            args.output,
        )
        print(f"Running command:\n\n{shlex.join(cmd)}", flush=True)
        subprocess.check_call(cmd)
        result = 0
    except (OSError, ValueError, subprocess.CalledProcessError) as e:
        print(e)
        result = 1
    if result != 0:
        print(f"An error occurded. Exiting with code {result}")
    return result


if __name__ == "__main__":
    sys.exit(main())
