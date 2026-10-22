"""Tests for java_options() in common.smk."""

import re
from pathlib import Path
from types import SimpleNamespace

COMMON = Path(__file__).parents[2] / "workflow/rules/common.smk"


def java_options(**kwargs):
    """Return the options that java_options(**kwargs) builds for a 4-thread job."""
    match = re.search(
        r"^def java_options\(.*?(?=^\S)", COMMON.read_text(), re.MULTILINE | re.DOTALL
    )
    assert match, "common.smk defines no java_options"
    scope = {"JAVA_TMP_DIR": ".picard"}
    exec(match.group(0), scope)  # noqa: S102, common.smk is not importable
    params = scope["java_options"](**kwargs)
    return params(wildcards=None, threads=4, resources=SimpleNamespace(mem_mb=10000))


def test_jvm_sizes_itself_to_the_job_threads():
    assert "-XX:ActiveProcessorCount=4" in java_options()


def test_heap_defaults_to_80_percent_of_mem_mb():
    assert "-Xmx8000m" in java_options()
    assert "-Xmx3000m" in java_options(heap_mb=3000)
