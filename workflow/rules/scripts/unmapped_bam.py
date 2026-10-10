import subprocess
import shlex
import gzip
import argparse


CMD = (
    """
    picard -Xms{mem_mb}m -Xmx{mem_mb}m -Dorg.xerial.snappy.tempdir={picard_tmp} -Djava.io.tmpdir={picard_tmp} FastqToSam \
    -FASTQ {fastq1} \
    -FASTQ2 {fastq2} \
    -OUTPUT {unmapped_bam} \
    -READ_GROUP_NAME {read_group_name} \
    -SAMPLE_NAME {sample_name} \
    -LIBRARY_NAME MOHD-{sample_name} \
    -PLATFORM_UNIT {platform_unit} \
    -PLATFORM {platform} \
    -SORT_ORDER queryname
    """
)

def parse_and_run(mem_mb, picard_tmp, sample, platform, read1, read2, unmapped_bam):
    try:
        with gzip.open(read1, "rt") as f:
            first_line = f.readline()
            whitespace = first_line.index(" ")
            [instrument_id, _, flowcell_id, lane_no, _, _, _, _] = first_line[:whitespace].lstrip("@").split(":")

        cmd = CMD.format(
            mem_mb=mem_mb,
            picard_tmp=picard_tmp,
            fastq1=read1,
            fastq2=read2,
            unmapped_bam=unmapped_bam,
            read_group_name=f"{flowcell_id}.{lane_no}",
            sample_name=sample,
            platform_unit=f"{instrument_id}",
            platform=platform
        )
        split_cmd = shlex.split(cmd)

        print(f"Running command:\n\n{cmd}")

        subprocess.check_call(split_cmd)
        return(0)
    except Exception as e:
        print(e)
        return(1)


def main(args):
    result = parse_and_run(args.mem_mb, args.picard_tmp, args.sample, args.platform, args.read1, args.read2, args.unmapped_bam)
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mem_mb", help="memory in mb")
    parser.add_argument("picard_tmp", help="picard tmp directory")
    parser.add_argument("sample", help="sample name")
    parser.add_argument("platform", help="platform")
    parser.add_argument("read1", help="read1 file")
    parser.add_argument("read2", help="read2 file")
    parser.add_argument("unmapped_bam", help="unmapped bam file")
    args = parser.parse_args()
    result = main(args)
    if result != 0:
        print(f"An error occurded. Exiting with code {result}")
    exit(result)
