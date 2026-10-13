"""Run bundled numerical scripts with bounded time and a clean work directory.

This is process isolation for trusted repository code, not a security sandbox.
New LLM-generated code execution is explicitly opt-in on a controlled host.
"""
from __future__ import annotations

import ast
from dataclasses import dataclass, field
from io import BytesIO
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Callable
import zipfile

from .archive import ROOT, redact

MAX_LOG = 160_000
EXECUTION_TIMEOUT_OPTIONS = (15, 30, 60, 120, 180, 300, 600, 1200, 1800)
DEFAULT_EXECUTION_TIMEOUT = 600
MAX_EXECUTION_TIMEOUT = max(EXECUTION_TIMEOUT_OPTIONS)


@dataclass
class ExecutionResult:
    stdout: str
    stderr: str
    returncode: int | None
    elapsed: float
    timed_out: bool
    images: dict[str, bytes] = field(default_factory=dict)
    bundle: bytes = b""

    @property
    def succeeded(self) -> bool:
        return not self.timed_out and self.returncode == 0


def live_enabled() -> bool:
    return os.getenv("RE4_ENABLE_LIVE_EXECUTION", "").strip().lower() in {"1", "true", "yes"}


def validate_code(code: str) -> None:
    """Catch inappropriate process/network operations before any launch.

    Deliberately not advertised as a complete untrusted-code sandbox.
    """
    tree = ast.parse(code)
    denied_imports = {"subprocess", "socket", "requests", "urllib", "http", "ftplib", "paramiko",
                      "ctypes", "multiprocessing", "shutil", "openai"}
    denied_calls = {"eval", "exec", "compile", "__import__", "system", "popen", "spawn", "fork",
                    "remove", "unlink", "rmdir", "rmtree", "getenv", "environ"}
    for node in ast.walk(tree):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            names = [n.name for n in node.names] if isinstance(node, ast.Import) else [node.module or ""]
            if any(name.split(".")[0] in denied_imports for name in names):
                raise ValueError("The execution adapter permits numerical scripts only; network/process imports are blocked.")
        if isinstance(node, ast.Call):
            name = node.func.id if isinstance(node.func, ast.Name) else node.func.attr if isinstance(node.func, ast.Attribute) else ""
            if name in denied_calls:
                raise ValueError(f"The execution adapter does not permit {name}().")
        if isinstance(node, ast.Attribute) and node.attr == "environ":
            raise ValueError("Reading environment variables is not permitted in generated scripts.")


def clean_environment(workspace: Path) -> dict[str, str]:
    # Do not forward credentials, proxy tokens, or cloud service secrets.
    allowed = {"PATH", "SYSTEMROOT", "WINDIR", "COMSPEC", "PATHEXT", "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH"}
    env = {k: v for k, v in os.environ.items() if k.upper() in allowed}
    bundled_runtime = (ROOT / ".runtime" / "python").resolve()
    library_paths = [p for p in sys.path if p and Path(p).is_dir()
                     and (not Path(p).resolve().is_relative_to(ROOT) or Path(p).resolve() == bundled_runtime)]
    env.update({"PYTHONPATH": os.pathsep.join(library_paths), "PYTHONIOENCODING": "utf-8",
                "PYTHONUNBUFFERED": "1", "MPLBACKEND": "Agg", "MPLCONFIGDIR": str(workspace / ".mpl"),
                "TMP": str(workspace),
                "TEMP": str(workspace), "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1"})
    return env


def execute_script(code: str, timeout: int = DEFAULT_EXECUTION_TIMEOUT, *, generated: bool = False,
                   on_progress: Callable[[str, str], None] | None = None) -> ExecutionResult:
    if generated and not live_enabled():
        raise PermissionError("Enable RE4_ENABLE_LIVE_EXECUTION on a controlled host to run newly generated code.")
    validate_code(code)
    timeout = max(5, min(int(timeout), MAX_EXECUTION_TIMEOUT))
    started = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="re4-lab-") as folder:
        work = Path(folder)
        (work / "solver.py").write_text(code, encoding="utf-8")
        for name in ("dataset_keyhole.csv", "dataset_rayben.csv"):
            source = ROOT / "src" / "pycodes_ver" / "workspace" / name
            if source.is_file():
                shutil.copy2(source, work / name)
        # Scientific limits are best-effort on Unix; the wall-clock limit also applies on Windows.
        launcher = f"""import os, runpy
try:
    import resource
    resource.setrlimit(resource.RLIMIT_CPU, ({timeout}, {timeout + 1}))
    resource.setrlimit(resource.RLIMIT_FSIZE, (64*1024*1024, 64*1024*1024))
except (ImportError, ValueError, OSError):
    pass
""" + """import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
def capture_open_figures(*args, **kwargs):
    for number in plt.get_fignums():
        name = f'captured_figure_{number}.png'
        plt.figure(number).savefig(name, dpi=150, bbox_inches='tight')
        print(f'[Re4 display adapter] Saved open Matplotlib figure: {name}', flush=True)
plt.show = capture_open_figures
try:
    runpy.run_path('solver.py', run_name='__main__')
finally:
    if not any(Path('.').glob('*.png')) and plt.get_fignums():
        capture_open_figures()
    plt.close('all')
"""
        (work / "_launch.py").write_text(launcher, encoding="utf-8")
        timed_out = False
        # Files bound host-side memory even when a script is excessively verbose.
        with (work / "stdout.txt").open("wb") as stdout_file, (work / "stderr.txt").open("wb") as stderr_file:
            process = subprocess.Popen([sys.executable, "-u", "_launch.py"], cwd=work,
                                       env=clean_environment(work), stdout=stdout_file, stderr=stderr_file)
            process_started = time.perf_counter()
            try:
                while process.poll() is None:
                    if on_progress:
                        stdout_file.flush()
                        stderr_file.flush()
                        preview_out = (work / "stdout.txt").read_bytes()[-MAX_LOG:].decode("utf-8", errors="replace")
                        preview_err = (work / "stderr.txt").read_bytes()[-MAX_LOG:].decode("utf-8", errors="replace")
                        on_progress(redact(preview_out), redact(preview_err))
                    if time.perf_counter() - process_started > timeout:
                        timed_out = True
                        break
                    if any((work / name).stat().st_size > 8_000_000 for name in ("stdout.txt", "stderr.txt")):
                        break
                    time.sleep(.35)
            finally:
                # Also stop a child when the UI session is interrupted mid-execution.
                if process.poll() is None:
                    process.kill()
                process.wait(timeout=10)
        def log(name: str) -> str:
            with (work / name).open("rb") as handle:
                data = handle.read(MAX_LOG)
            text = data.decode("utf-8", errors="replace")
            if (work / name).stat().st_size > MAX_LOG:
                text += "\n[Display truncated at 160 kB]"
            return redact(text)
        stdout, stderr = log("stdout.txt"), log("stderr.txt")
        images: dict[str, bytes] = {}
        downloadable: dict[str, bytes] = {"solver.py": code.encode(), "stdout.txt": stdout.encode(), "stderr.txt": stderr.encode()}
        for path in sorted(work.rglob("*")):
            if not path.is_file() or path.is_symlink() or ".mpl" in path.parts:
                continue
            name = path.relative_to(work).as_posix()
            if path.suffix.lower() in {".png", ".jpg", ".jpeg", ".gif", ".webp"} and path.stat().st_size < 16_000_000:
                images[name] = path.read_bytes()
                downloadable[name] = images[name]
            elif path.suffix.lower() in {".csv", ".txt", ".json", ".pdf"} and path.stat().st_size < 2_000_000 and name not in {"stdout.txt", "stderr.txt"}:
                downloadable[name] = path.read_bytes()
        buffer = BytesIO()
        with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as bundle:
            for name, data in downloadable.items():
                bundle.writestr(name, data)
        return ExecutionResult(stdout, stderr, process.returncode, time.perf_counter() - started,
                               timed_out, images, buffer.getvalue())
