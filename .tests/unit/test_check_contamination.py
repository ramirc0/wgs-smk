"""Tests for check_contamination.py."""

import check_contamination
import pytest

HEADER = "#SEQ_ID\tRG\tCHIP_ID\t#SNPS\t#READS\tAVG_DP\tFREEMIX\tFREELK1\tFREELK0\n"
ROW = "S\tNA\tNA\t100000\t2823125\t28.2\t0.0015\t-591218\t-593292\n"


@pytest.mark.parametrize("rows, expected", [(1, 0), (0, 1), (2, 1)])
def test_check_self_sm_requires_one_row(tmp_path, rows, expected):
    self_sm = tmp_path / "S.selfSM"
    self_sm.write_text(HEADER + ROW * rows)
    assert check_contamination.check_self_sm(self_sm, 0.75) == expected


def test_check_self_sm_rejects_zero_likelihoods(tmp_path):
    self_sm = tmp_path / "S.selfSM"
    self_sm.write_text(HEADER + "S\tNA\tNA\t0\t0\t0\t0\t0\t0\n")
    assert check_contamination.check_self_sm(self_sm, 0.75) == 1


@pytest.mark.parametrize("name", ["S.preBqsr", "Ms"])
def test_verifybamid_output_prefix_drops_only_the_suffix(tmp_path, monkeypatch, name):
    calls = []
    monkeypatch.setattr(check_contamination.subprocess, "check_call", calls.append)
    output = tmp_path / f"{name}.selfSM"
    output.write_text(HEADER + ROW)
    argv = ["--threads", "4", "--output", str(output), "--bam", "b", "--reference"]
    argv += ["r", "--svd-prefix", "p", "--underestimation-factor"]
    assert check_contamination.main([*argv, "0.75"]) == 0
    cmd = calls[0]
    assert cmd[cmd.index("--Output") + 1] == str(tmp_path / name)
