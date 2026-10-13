"""Tests for unmapped_bam.py read group routing."""

import gzip
import subprocess

import pytest
import unmapped_bam

READS = [
    ("NS500269:341:HKGGNBGXY:1:11101:13293:1046", "1"),
    ("LH00180:355:23CWFWLT4:3:1101:4964:1000:CCGTAATTG", "GTACTCGTGA+CACTTGAGAA"),
    ("NS500269:341:HKGGNBGXY:2:11101:6601:1046", "1"),
    ("LH00180:355:23CWFWLT4:3:1101:5000:1000:TTGCCAATG", "GTACTCGTGA+CACTTGAGAA"),
]


def write_pair(tmp_path, reads):
    fastqs = []
    for end in (1, 2):
        path = tmp_path / f"R{end}.fastq.gz"
        with gzip.open(path, "wt") as f:
            for name, barcode in reads:
                f.write(f"@{name} {end}:N:0:{barcode}\nACGT\n+\nIIII\n")
        fastqs.append(str(path))
    return fastqs


def test_route_splits_7_and_8_field_names_by_flowcell_and_lane(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    read_groups = unmapped_bam.route(*write_pair(tmp_path, READS), str(out))
    assert read_groups == {
        "HKGGNBGXY.1": ("1", 1),
        "HKGGNBGXY.2": ("1", 1),
        "23CWFWLT4.3": ("GTACTCGTGA+CACTTGAGAA", 2),
    }
    r2 = (out / "23CWFWLT4.3_R2.fastq").read_text().splitlines()
    assert r2[::4] == [f"@{READS[i][0]} 2:N:0:{READS[i][1]}" for i in (1, 3)]


def test_route_takes_the_most_common_barcode_of_a_read_group(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    reads = [
        ("LH00180:355:23CVNHLT4:8:1101:1000:1000:CCGTAATTG", "GNACTCGTGA+CACTTGAGAA"),
        ("LH00180:355:23CVNHLT4:8:1101:2000:1000:TTGCCAATG", "GTACTCGTGA+CACTTGAGAA"),
        ("LH00180:355:23CVNHLT4:8:1101:3000:1000:AACCGGTTA", "GTACTCGTGA+CACTTGAGAA"),
    ]
    read_groups = unmapped_bam.route(*write_pair(tmp_path, reads), str(out))
    assert read_groups == {"23CVNHLT4.8": ("GTACTCGTGA+CACTTGAGAA", 3)}


def test_route_rejects_other_name_shapes(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    with pytest.raises(subprocess.CalledProcessError):
        unmapped_bam.route(*write_pair(tmp_path, [("A:B:C:D:E:F", "1")]), str(out))
