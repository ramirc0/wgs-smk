"""Convert paired FASTQs to one unmapped BAM per read group, as WARP takes them.

Each read pair goes to the read group of its flowcell and lane, read name fields
3 and 4. FastqToSam converts the read groups in parallel. A read group uBAM
larger than `--cutoff-gib` is also split into files of `--reads-per-file` reads,
as WARP's SplitLargeReadGroup does.
"""

import argparse
import os
import shlex
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

# Input: one tab-joined R1 record and R2 record per line. Writes each pair to
# <dir>/<flowcell>.<lane>_R{1,2}.fastq. The last field of the first R1
# comment of a read group is its barcode. Prints one line per read group:
# ID, barcode, pairs.
ROUTE = r"""
{
    words = split($1, word, " ")
    fields = split(word[1], field, ":")
    if (fields != 7 && fields != 8) {
        print "Expected 7 or 8 colon-separated fields: " $1 > "/dev/stderr"
        exit 1
    }
    rg = field[3] "." field[4]
    if (!(rg in pairs)) {
        n = words > 1 ? split(word[2], comment, ":") : 0
        barcode[rg] = n ? comment[n] : ""
    }
    pairs[rg]++
    print $1 "\n" $2 "\n" $3 "\n" $4 > (dir "/" rg "_R1.fastq")
    print $5 "\n" $6 "\n" $7 "\n" $8 > (dir "/" rg "_R2.fastq")
}
END { for (rg in pairs) print rg "\t" barcode[rg] "\t" pairs[rg] }
"""


def build_parser():
    """Build the unmapped_bam.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the unmapped_bam.py flags.
    """
    p = argparse.ArgumentParser(description="Convert paired FASTQs to uBAMs.")
    p.add_argument("--java-options", required=True, help="Options of each JVM.")
    p.add_argument("--tmp-dir", required=True, help="Scratch for FASTQs and spills.")
    p.add_argument("--jobs", type=int, required=True, help="Concurrent Picard runs.")
    p.add_argument("--sample", required=True, help="Sample ID.")
    p.add_argument("--platform", required=True, help="Read group platform.")
    p.add_argument("--r1", required=True, help="Gzipped R1 FASTQ.")
    p.add_argument("--r2", required=True, help="Gzipped R2 FASTQ.")
    p.add_argument(
        "--cutoff-gib", type=float, required=True, help="Split uBAMs above this size."
    )
    p.add_argument(
        "--reads-per-file", type=int, required=True, help="Reads per split uBAM."
    )
    p.add_argument(
        "--compression-level", type=int, required=True, help="Shard compression."
    )
    p.add_argument("--outdir", required=True, help="Directory for the uBAMs.")
    return p


def route(r1, r2, fastq_dir):
    """Split a FASTQ pair into one FASTQ pair per read group.

    Parameters
    ----------
    r1, r2 : str
        Gzipped R1 and R2 FASTQs. Every read name MUST hold 7 (Casava 1.8) or
        8 (with a UMI) colon-separated fields.
    fastq_dir : str
        Directory for the `<read group>_R1.fastq` and `_R2.fastq` files.

    Returns
    -------
    dict of str to tuple of (str, int)
        Read group ID, `<flowcell>.<lane>`, to its barcode and read pairs.
    """
    # Process substitution would hide a failed gzip, so use FIFOs and wait.
    fifos = [Path(fastq_dir, f"R{end}.fifo") for end in (1, 2)]
    readers = []
    for fastq, fifo in zip((r1, r2), fifos):
        os.mkfifo(fifo)
        readers.append(
            subprocess.Popen(
                ["bash", "-o", "pipefail", "-c", 'gzip -dc "$0" | paste - - - - > "$1"']
                + [fastq, str(fifo)]
            )
        )
    out = subprocess.run(
        [
            "bash",
            "-o",
            "pipefail",
            "-c",
            'paste "$0" "$1" | mawk -F "\t" -v dir="$2" "$3"',
        ]
        + [*map(str, fifos), fastq_dir, ROUTE],
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    ).stdout
    for reader in readers:
        if reader.wait():
            raise subprocess.CalledProcessError(reader.returncode, reader.args)
    for fifo in fifos:
        fifo.unlink()
    read_groups = {}
    for line in out.splitlines():
        rg, barcode, pairs = line.split("\t")
        read_groups[rg] = (barcode, int(pairs))
    return read_groups


def picard(java_options, tool, *args):
    """Run a Picard tool and raise if it fails.

    Parameters
    ----------
    java_options : str
        JVM options, split like a shell command line.
    tool : str
        Picard tool name.
    *args : str
        Tool arguments.
    """
    cmd = ["picard", *shlex.split(java_options), tool, *args]
    print(f"Running command:\n\n{shlex.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True)


def convert(args, fastq_dir, rg, barcode, pairs):
    """Write one read group's uBAM, and its shards if it is large.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command line.
    fastq_dir : str
        Directory with the read group's FASTQ pair. Both files are deleted.
    rg : str
        Read group ID, `<flowcell>.<lane>`.
    barcode : str
        Read group barcode. Empty if the FASTQ names have no comment.
    pairs : int
        Read pairs in the read group.
    """
    fastqs = [Path(fastq_dir, f"{rg}_R{end}.fastq") for end in (1, 2)]
    ubam = Path(args.outdir, f"{rg}.bam")
    picard(
        args.java_options,
        "FastqToSam",
        "-FASTQ",
        str(fastqs[0]),
        "-FASTQ2",
        str(fastqs[1]),
        "-OUTPUT",
        str(ubam),
        "-READ_GROUP_NAME",
        rg,
        "-SAMPLE_NAME",
        args.sample,
        "-LIBRARY_NAME",
        args.sample,
        "-PLATFORM_UNIT",
        f"{rg}.{barcode}" if barcode else rg,
        "-PLATFORM",
        args.platform,
        "-SORT_ORDER",
        "queryname",
        "-TMP_DIR",
        args.tmp_dir,
        # 10x the default: fewer, larger spills.
        "-MAX_RECORDS_IN_RAM",
        "5000000",
    )
    for fastq in fastqs:
        fastq.unlink()
    if ubam.stat().st_size > args.cutoff_gib * 2**30:
        shards = Path(args.outdir, rg)
        shards.mkdir()
        picard(
            f"{args.java_options} -Dsamjdk.compression_level={args.compression_level}",
            "SplitSamByNumberOfReads",
            "-INPUT",
            str(ubam),
            "-OUTPUT",
            str(shards),
            "-SPLIT_TO_N_READS",
            str(args.reads_per_file),
            "-TOTAL_READS_IN_INPUT",
            str(2 * pairs),
        )


def main(argv=None):
    """Route the reads by read group, then convert each read group.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    Path(args.outdir).mkdir(parents=True)
    with tempfile.TemporaryDirectory(dir=args.tmp_dir) as fastq_dir:
        read_groups = route(args.r1, args.r2, fastq_dir)
        print(f"Read groups: {read_groups}", flush=True)
        with ThreadPoolExecutor(args.jobs) as pool:
            for future in [
                pool.submit(convert, args, fastq_dir, rg, barcode, pairs)
                for rg, (barcode, pairs) in sorted(read_groups.items())
            ]:
                future.result()


if __name__ == "__main__":
    sys.exit(main())
