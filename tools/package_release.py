"""Build a clean uploadable source package, excluding local runtimes and credentials."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import zipfile

ROOT = Path(__file__).resolve().parents[1]
EXCLUDED = {'.git', '.venv', 'venv', '.runtime', '.showcase_work', '_showcase_work',
            '__pycache__', '.pytest_cache', '.ipynb_checkpoints', '.mypy_cache', 'runs', 'tmp'}


def files():
    for directory, directories, names in os.walk(ROOT):
        directories[:] = sorted(x for x in directories if x not in EXCLUDED)
        for name in sorted(names):
            path = Path(directory) / name
            if not path.is_file() or path.is_symlink():
                continue
            if name in {'secrets.toml', '.DS_Store', 'Thumbs.db'}:
                continue
            if name.startswith('.env') and name != '.env.example':
                continue
            if path.suffix.lower() in {'.pyc', '.pyo', '.log', '.zip'}:
                continue
            yield path


def main(destination: Path):
    chosen = list(files())
    leaks = []
    for path in chosen:
        if path.suffix in {'.py', '.md', '.ipynb', '.json', '.toml', '.txt', '.yaml', '.yml'}:
            text = path.read_text(encoding='utf-8-sig', errors='replace')
            # Ignore tests' deliberately synthetic strings; credentials in source are never allowed.
            if path.relative_to(ROOT).parts[0] != 'tests' and re.search(r'\b(?:sk-[A-Za-z0-9_-]{20,}|ghp_[A-Za-z0-9]{30,}|hf_[A-Za-z0-9]{30,})', text):
                leaks.append(path.relative_to(ROOT).as_posix())
    if leaks:
        raise SystemExit('Credential-shaped values found in: ' + ', '.join(leaks))
    destination.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(destination, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for path in chosen:
            archive.write(path, Path('RE4_code') / path.relative_to(ROOT))
    with zipfile.ZipFile(destination) as archive:
        if archive.testzip() is not None:
            raise SystemExit('Archive verification failed.')
    summary = {'archive': destination.name, 'files': len(chosen), 'bytes': destination.stat().st_size,
               'sha256': hashlib.sha256(destination.read_bytes()).hexdigest(),
               'excluded': sorted(EXCLUDED), 'credential_scan': 'passed'}
    destination.with_suffix('.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT.parent / 'output' / 'RE4_Streamlit_GitHub_ready.zip')
    main(parser.parse_args().output.resolve())
