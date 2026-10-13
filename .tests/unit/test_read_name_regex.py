"""Tests for MarkDuplicates' READ_NAME_REGEX in align.smk."""

import re
from pathlib import Path

import pytest

ALIGN = Path(__file__).parents[2] / "workflow/rules/align.smk"


def read_name_regex():
    """Return the READ_NAME_REGEX that mark_duplicates passes to Picard."""
    match = re.search(r'read_name_regex=lambda w: r"([^"]+)"', ALIGN.read_text())
    assert match, "mark_duplicates sets no read_name_regex"
    return match.group(1)


@pytest.mark.parametrize(
    "name",
    [
        "NS500269:341:HKGGNBGXY:1:1101:1558:1014",
        "LH00195:225:22MHVFLT4:2:1101:1558:1014:GTGAACGAA",
    ],
)
def test_regex_takes_tile_x_y_from_7_and_8_field_names(name):
    # Picard's ReadNameParser calls matches(), a full match.
    match = re.fullmatch(read_name_regex(), name)
    assert match and match.groups() == ("1101", "1558", "1014")


def test_regex_rejects_other_shapes():
    assert not re.fullmatch(read_name_regex(), "A:B:1101:1558:1014:X:Y:Z:W")
