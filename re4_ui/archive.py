"""Read actual archived agent reports without importing the LLM configuration."""
from __future__ import annotations

from dataclasses import dataclass
import ast
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE_ROOT = ROOT / "src" / "ipynb_ver" / "sample_results"
LABELS = {
    "sodShock_problem": "Sod shock tube",
    "cavityFlow_problem": "Lid-driven cavity",
    "cylinder_problem": "Cylinder flow",
    "longTmNS_problem": "Unsteady Navier–Stokes",
    "burgers_problem": "Burgers equation",
    "poison_problem": "Poisson equation",
    "helmholtz_problem": "Helmholtz equation",
}


def redact(text: str) -> str:
    """Prevent historical credentials from appearing in reports or downloads."""
    text = re.sub(r"\bsk-[A-Za-z0-9_-]{12,}", "[REDACTED]", text)
    return re.sub(
        r"(?im)^([^\n]*(?:api[_ -]?key|authorization|bearer)[^\n]*?[:=])[^\n]*$",
        r"\1 [REDACTED]", text,
    )


def sections(text: str) -> dict[str, str]:
    """Split level-three report headings, ignoring headings inside code fences."""
    result: dict[str, list[str]] = {}
    key = "preamble"
    result[key] = []
    fenced = False
    for line in text.splitlines():
        if line.strip().startswith("```"):
            fenced = not fenced
        known = {"prob_todo", "expanded_prob", "solution_plans", "technical_spec", "python_codes",
                 "runtime_outputs", "review_decision", "review_comments", "iteration_count", "rev_image_description"}
        if not fenced and line.startswith("### ") and line[4:].strip() in known:
            key = line[4:].strip()
            result.setdefault(key, [])
        else:
            result[key].append(line)
    return {k: "\n".join(v).strip() for k, v in result.items()}


def stage_blocks(text: str) -> list[str]:
    if not text.strip():
        return []
    items = re.split(r"(?:\s*Current Stage\s*\[[A-Z]/\d+\]\s*|(?m:^#### (?:Output|Script) block\s*\d+:?[^\n]*\n))", text)
    return [item.strip() for item in items if item.strip()]


def python_blocks(text: str) -> list[str]:
    return [m.strip() for m in re.findall(r"```(?:python|py)\s*\n(.*?)```", text, re.S | re.I)]


def report_codes(text: str) -> list[str]:
    """Some historical reports saved bare Python instead of Markdown fences."""
    blocks = re.split(r"(?m)^#### Script block\s*\d+:?[^\n]*\n", text)
    if len(blocks) == 1:
        return python_blocks(text)
    result = []
    for block in blocks[1:]:
        fenced = python_blocks(block)
        result.append(fenced[0] if fenced else block.strip())
    return result


def solution_plans(text: str) -> list[tuple[str, str]]:
    """Decode the original Pydantic string representation for readable cards."""
    literal = r"(?:'(?:\\.|[^'\\])*'|\"(?:\\.|[^\"\\])*\")"
    pattern = rf"solu_name=({literal})\s+content=({literal})"
    plans = []
    for title, content in re.findall(pattern, text, re.S):
        try:
            plans.append((ast.literal_eval(title), ast.literal_eval(content)))
        except (ValueError, SyntaxError):
            continue
    return plans


@dataclass(frozen=True)
class Round:
    index: int
    code: str
    specification: str = ""
    runtime: str = ""
    review: str = ""
    image_review: str = ""


@dataclass(frozen=True)
class Archive:
    path: Path
    title: str
    problem: str
    expanded_problem: str
    plans: str
    rounds: tuple[Round, ...]
    final_decision: str
    raw: str


def archive_choices() -> dict[str, Path]:
    choices: dict[str, Path] = {}
    for path in sorted(ARCHIVE_ROOT.rglob("*.md")):
        if path.parent.name not in LABELS:
            continue
        name = path.stem.split("_wt_")[-1]
        choices[f"{LABELS[path.parent.name]} · trial {name}"] = path
    return choices


def read_archive(path: str | Path) -> Archive:
    path = Path(path).resolve()
    if not path.is_relative_to(ARCHIVE_ROOT.resolve()) or path.suffix.lower() != ".md":
        raise ValueError("Choose a Markdown report from the bundled sample_results archive.")
    raw = redact(path.read_text(encoding="utf-8-sig", errors="replace"))
    parts = sections(raw)
    codes = report_codes(parts.get("python_codes", ""))
    specs = stage_blocks(parts.get("technical_spec", ""))
    outputs = stage_blocks(parts.get("runtime_outputs", ""))
    reviews = stage_blocks(parts.get("review_comments", ""))
    images = stage_blocks(parts.get("rev_image_description", ""))
    def at(items: list[str], i: int) -> str:
        return items[i] if i < len(items) else ""
    rounds = tuple(Round(i, code, at(specs, i), at(outputs, i), at(reviews, i), at(images, i))
                   for i, code in enumerate(codes))
    return Archive(path, LABELS.get(path.parent.name, path.stem),
                   parts.get("prob_todo", ""), parts.get("expanded_prob", ""),
                   parts.get("solution_plans", ""), rounds,
                   parts.get("review_decision", "not recorded").strip(), raw)


def code_choices() -> dict[str, tuple[Path, int]]:
    choices: dict[str, tuple[Path, int]] = {}
    for path in sorted(ARCHIVE_ROOT.rglob("solu_rev_2.dat")):
        relative = path.relative_to(ARCHIVE_ROOT)
        choices[" / ".join(relative.parts[:-1])] = (path, 0)
    return choices


def read_trusted_code(path: Path, index: int = 0) -> str:
    path = path.resolve()
    if not path.is_relative_to(ARCHIVE_ROOT.resolve()) or path.suffix.lower() not in {".md", ".dat"}:
        raise ValueError("Only bundled archived scripts can be selected.")
    if path.suffix.lower() == ".md":
        return read_archive(path).rounds[index].code
    blocks = python_blocks(path.read_text(encoding="utf-8-sig", errors="replace"))
    if not blocks:
        raise ValueError("This archive contains no fenced Python script.")
    return redact(blocks[index])
