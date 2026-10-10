import argparse
import subprocess
import glob 
import os
import sys

CMD = (
    """
    picard -Xms{mem_mb}m -Xmx{mem_mb}m \
    IntervalListTools \
          -SCATTER_COUNT {scatter_count} \
          -SUBDIVISION_MODE BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
          -UNIQUE true \
          -SORT true \
          -BREAK_BANDS_AT_MULTIPLES_OF {break_bands_at_multiples_of} \
          -INPUT {interval_list} \
          -OUTPUT {outdir}
    """
)

def main(args):
    return_code = 0

    try:
        mem_mb = args.mem_mb
        scatter_count = args.scatter_count
        break_bands_at_multiples_of = args.break_bands_at_multiples_of
        interval_list = args.interval_list
        outdir = args.outdir

        picard_command = CMD.format(
            mem_mb=mem_mb,
            scatter_count=scatter_count,
            break_bands_at_multiples_of=break_bands_at_multiples_of,
            interval_list=interval_list,
            outdir=outdir
        )

        subprocess.check_call(picard_command, shell=True)
    except Exception as e:
        return_code = 1
        print(e)
        return return_code
    finally:
        outdir = args.outdir

        intervals = sorted(glob.glob(f"{outdir}/*/*.interval_list"))
        for i, interval in enumerate(intervals):
            (directory, filename) = os.path.split(interval)
            newName = os.path.join(directory, str(i + 1) + filename)
            os.rename(interval, newName)
        print(len(intervals))
        return return_code

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mem_mb", help="memory in mb")
    parser.add_argument("scatter_count", type=int, help="scatter count")
    parser.add_argument("break_bands_at_multiples_of", type=int, help="break bands at multiples of")
    parser.add_argument("interval_list", help="interval list")
    parser.add_argument("outdir", help="output directory")

    args = parser.parse_args()
    result = main(args)
    if result != 0:
        print(f"An error occurded. Exiting with code {result}")
    sys.exit(result)

