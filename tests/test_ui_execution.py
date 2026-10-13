from io import BytesIO
from pathlib import Path
import zipfile

import pytest

from re4_ui.archive import ARCHIVE_ROOT, read_trusted_code
from re4_ui.execution import clean_environment, execute_script, validate_code


def test_numerical_archive_executes_and_exports_its_open_figure():
    path = next(p for p in ARCHIVE_ROOT.rglob("solu_rev_2.dat")
                if "SodShock" in p.parts and "GPT_ans" in p.parts and "sample1" in p.parts)
    progress = []
    result = execute_script(read_trusted_code(path), timeout=30,
                            on_progress=lambda out, err: progress.append((out, err)))
    assert result.succeeded, result.stderr
    assert "Solution summary" in result.stdout
    assert len(result.images) >= 1
    assert progress
    with zipfile.ZipFile(BytesIO(result.bundle)) as archive:
        assert "solver.py" in archive.namelist()
        assert "stdout.txt" in archive.namelist()
        assert any(n.endswith(".png") for n in archive.namelist())


def test_runtime_failure_is_preserved():
    result = execute_script("print('before failure', flush=True)\nraise RuntimeError('numerical failure')", timeout=15)
    assert not result.succeeded
    assert result.returncode != 0
    assert "before failure" in result.stdout
    assert "RuntimeError: numerical failure" in result.stderr


def test_wall_clock_timeout_is_enforced():
    result = execute_script("while True: pass", timeout=5)
    assert result.timed_out
    assert not result.succeeded
    assert result.elapsed < 12


def test_child_environment_omits_keys_and_retains_dependencies(monkeypatch, tmp_path):
    monkeypatch.setenv("OPENAI_API_KEY", "example-private-value")
    monkeypatch.setenv("RE4_API_KEY", "another-private-value")
    env = clean_environment(tmp_path)
    assert "OPENAI_API_KEY" not in env
    assert "RE4_API_KEY" not in env
    assert "example-private-value" not in repr(env)
    assert env["MPLCONFIGDIR"] == str(tmp_path / ".mpl")


@pytest.mark.parametrize("code", ["import requests", "import subprocess", "import os\nos.system('whoami')", "eval('1+1')", "import os\nprint(os.environ)"])
def test_process_network_and_environment_operations_are_rejected(code):
    with pytest.raises(ValueError):
        validate_code(code)


def test_generated_execution_is_disabled_by_default(monkeypatch):
    monkeypatch.delenv("RE4_ENABLE_LIVE_EXECUTION", raising=False)
    with pytest.raises(PermissionError):
        execute_script("print('not run')", generated=True)
