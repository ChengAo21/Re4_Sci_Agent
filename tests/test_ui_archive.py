from pathlib import Path
import re

import pytest

from re4_ui.archive import ARCHIVE_ROOT, archive_choices, read_archive, read_trusted_code, redact, sections, solution_plans


def test_all_archived_reports_preserve_every_attempt():
    choices = archive_choices()
    assert len(choices) == 20
    for path in choices.values():
        raw = path.read_text(encoding="utf-8-sig")
        report = read_archive(path)
        count = len(re.findall(r"(?m)^#### Script block\d+", raw))
        assert len(report.rounds) == count, str(path)
        assert report.problem
        if "ConsModule_ON: False" not in raw:
            assert report.expanded_problem and report.plans
            assert len(solution_plans(report.plans)) == report.plans.count("solu_name=")
        assert all(r.code and r.runtime and r.review for r in report.rounds), str(path)


def test_parser_does_not_split_python_headings():
    source = "### python_codes\n```python\n### runtime_outputs\nprint('a')\n```\n### runtime_outputs\nreal output"
    result = sections(source)
    assert "print('a')" in result["python_codes"]
    assert result["runtime_outputs"] == "real output"


def test_paths_outside_bundled_archive_are_rejected(tmp_path):
    with pytest.raises(ValueError):
        read_archive(tmp_path / "report.md")
    with pytest.raises(ValueError):
        read_trusted_code(tmp_path / "script.dat")


def test_credentials_are_redacted():
    sample = "API_KEY=example-private-value\nAuthorization: Bearer abcdefg\nkey: sk-1234567890abcdefghijklmnop\nordinary scientific output"
    result = redact(sample)
    assert "example-private-value" not in result
    assert "abcdefg" not in result
    assert "sk-123" not in result
    assert "ordinary scientific output" in result
