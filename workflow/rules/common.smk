# Config, sample sheet, and shared helpers.

from pathlib import Path

import polars as pl
from snakemake.exceptions import WorkflowError


configfile: "config/config.yaml"


RESULTS = config["outdir"]
RUN_ID = config.get("run_id", "default")
OUTDIR = f"{RESULTS}/{RUN_ID}"          # per-run
REFDIR = f"{RESULTS}/refs"              # run-independent
LOGDIR = f"logs/{RUN_ID}"
BENCHDIR = f"benchmarks/{RUN_ID}"
REFS = config["references"]

# The index files sit beside the FASTA. Every rule that reads the FASTA lists
# all three, so reference_index runs first.
REFERENCE_FASTA = REFS["fasta"]
REFERENCE = {
    "fasta": REFERENCE_FASTA,
    "fasta_index": f"{REFERENCE_FASTA}.fai",
    "fasta_dict": str(Path(REFERENCE_FASTA).with_suffix(".dict")),
}
CHROMSIZES = str(Path(REFERENCE_FASTA).with_suffix(".chromsizes"))
DRAGEN_REFERENCE = f"{REFDIR}/dragen_reference"

SCATTER_COUNT = config["scatter_interval_list"]["scatter_count"]
SCATTERS = range(1, SCATTER_COUNT + 1)
SCATTERED_INTERVAL_LIST = (
    f"{OUTDIR}/scatter_interval_list/temp_{{scatter:04d}}_of_{SCATTER_COUNT}"
    "/{scatter}scattered.interval_list"
)

# The workspace path is relative to the rule file. Every rule file sits in rules/.
SOFTWARE_ENV = pixi(workspace="../envs", env=config["pixi_env"], locked=True)

# GKL and snappy extract native libraries to java.io.tmpdir, and a noexec /tmp
# blocks them. Picard moves java.io.tmpdir to its TMP_DIR only after GKL loads.
JAVA_TMP_DIR = config.get("picard_tmp_dir", ".picard")
Path(JAVA_TMP_DIR).mkdir(parents=True, exist_ok=True)


_SHEET = Path(config["samples"])
if not _SHEET.exists():
    raise WorkflowError(
        f"Sample sheet not found: {_SHEET}\n"
        "Columns: sample_id, r1, r2; see config/samples.example.tsv."
    )

# All-string read so empty cells are "" (absent), not a type error.
_manifest = pl.read_csv(_SHEET, separator="\t", infer_schema_length=0, comment_prefix="#")


def _column(name):
    """Map each sample to its value in sheet column `name`.

    Parameters
    ----------
    name : str
        Sample sheet column.

    Returns
    -------
    dict of str to str or None
        Value per `sample_id`, with empty cells as None. Empty when the sheet
        has no such column.
    """
    if name not in _manifest.columns:
        return {}
    return {
        row["sample_id"]: (row[name] or None)
        for row in _manifest.iter_rows(named=True)
    }


SAMPLES = _manifest.get_column("sample_id").to_list()
R1_OF = _column("r1")
R2_OF = _column("r2")

if not SAMPLES:
    raise WorkflowError(f"Sample sheet {_SHEET} has no rows.")


# One BQSR group per non-empty line: the tab-separated contigs it covers.
SEQUENCE_GROUPS = [
    line.strip().split("\t")
    for line in Path(REFS["sequence_grouping"]).read_text().strip().split("\n")
    if line.strip()
]
GROUPS = range(len(SEQUENCE_GROUPS))


wildcard_constraints:
    sample=r"[A-Za-z0-9][A-Za-z0-9_.-]*",
    group=r"\d+",
    scatter=r"\d+",


def repeat_flag(flag, values):
    """Return `flag` before each value, as one token list.

    Parameters
    ----------
    flag : str
        Command-line flag, such as `-L` or `-INPUT`.
    values : iterable
        Values to pass. Each is converted with `str`.

    Returns
    -------
    list of str
        `[flag, value1, flag, value2, ...]`. Interpolate it with `:q`.
    """
    return [token for value in values for token in (flag, str(value))]


def java_options(heap_mb=None):
    """Return a `params:` function that builds a Java rule's JVM options.

    Parameters
    ----------
    heap_mb : int, optional
        Maximum heap in MB. Defaults to 80% of the job's `mem_mb`.

    Returns
    -------
    callable
        Maps the job's wildcards and resources to the options, as a list.
        Interpolate it with `:q`, or bare inside a double-quoted string.
    """

    def options(wildcards, resources):
        heap = heap_mb or int(resources.mem_mb * 0.8)
        return [
            f"-Xmx{heap}m",
            f"-Djava.io.tmpdir={JAVA_TMP_DIR}",
            f"-Dorg.xerial.snappy.tempdir={JAVA_TMP_DIR}",
        ]

    return options
