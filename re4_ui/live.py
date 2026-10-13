"""Stream the real Re4 stages through an explicitly configured API.

Original schema classes and static system prompts are read from src/pycodes_ver
without importing its notebook-oriented runner or credential configuration.
"""
from __future__ import annotations

import ast
from dataclasses import dataclass
import json
from pathlib import Path
import time
from typing import Any, Callable

from .archive import ROOT, python_blocks, redact
from .execution import DEFAULT_EXECUTION_TIMEOUT, execute_script, live_enabled, ExecutionResult


@dataclass
class LiveSettings:
    api_key: str
    base_url: str
    consultant_model: str
    programmer_model: str
    reviewer_model: str
    revisions: int = 2
    execution_timeout: int = DEFAULT_EXECUTION_TIMEOUT
    max_tokens: int = 10000
    image_review: bool = False


def problems() -> dict[str, str]:
    source = (ROOT / "src" / "pycodes_ver" / "problems.py").read_text(encoding="utf-8-sig")
    result = {}
    for node in ast.parse(source).body:
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    result[target.id] = node.value.value.strip()
    return result


def original_contract() -> tuple[dict[str, Any], dict[str, str]]:
    from pydantic import BaseModel, Field
    from typing import List, Literal
    schema_path = ROOT / "src" / "pycodes_ver" / "schemas.py"
    tree = ast.parse(schema_path.read_text(encoding="utf-8-sig"))
    definitions = [n for n in tree.body if isinstance(n, ast.ClassDef) and any(isinstance(b, ast.Name) and b.id == "BaseModel" for b in n.bases)]
    # Only the repository's Pydantic definitions are executed. No config module is imported.
    namespace = {"BaseModel": BaseModel, "Field": Field, "List": List, "Literal": Literal,
                 "ALTERNATIVE_NUM": 2, "MAX_COMMENTS": 5, "__name__": "re4_ui.live_contract"}
    exec(compile(ast.Module(body=definitions, type_ignores=[]), str(schema_path), "exec", dont_inherit=True), namespace)
    node_tree = ast.parse((ROOT / "src" / "pycodes_ver" / "nodes.py").read_text(encoding="utf-8-sig"))
    wanted = {"sys_msg", "sys_msg_coder", "add_first_sys_msg", "add_rev_sys_msg", "sys_msg_reviewer", "image_rev_sys_msg"}
    prompts = {}
    for node in ast.walk(node_tree):
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in wanted:
                    # add_first_sys_msg has a no-Consultant variant later: retain first one.
                    prompts.setdefault(target.id, node.value.value)
    return namespace, prompts


def run_pipeline(problem: str, settings: LiveSettings, emit: Callable[[dict], None]) -> list[dict]:
    if not live_enabled():
        raise PermissionError("Live execution is disabled on this host.")
    if not settings.api_key.strip() or not all([settings.consultant_model, settings.programmer_model, settings.reviewer_model]):
        raise ValueError("Provide an API key and model names for all three roles.")
    from openai import OpenAI
    schemas, prompts = original_contract()
    client = OpenAI(api_key=settings.api_key, base_url=settings.base_url or None, timeout=120.0, max_retries=0)
    events: list[dict] = []

    def publish(event: dict) -> None:
        events.append(event)
        emit(event)

    def call(role: str, model: str, system: str, user: Any, schema_name: str, iteration: int) -> dict:
        schema = schemas[schema_name]
        system += "\nReturn a JSON object matching this JSON schema exactly:\n" + json.dumps(schema.model_json_schema())
        started = time.perf_counter()
        publish({"type": "start", "role": role, "iteration": iteration, "model": model})
        pieces = []
        usage = None
        stream = client.chat.completions.create(
            model=model, messages=[{"role": "system", "content": system}, {"role": "user", "content": user}],
            response_format={"type": "json_object"}, stream=True, stream_options={"include_usage": True},
            max_completion_tokens=settings.max_tokens,
        )
        for chunk in stream:
            if chunk.usage:
                usage = chunk.usage.model_dump()
            if chunk.choices and chunk.choices[0].delta.content:
                pieces.append(chunk.choices[0].delta.content)
                emit({"type": "token", "role": role, "iteration": iteration, "text": redact("".join(pieces))})
        content = "".join(pieces)
        response = schema.model_validate_json(content).model_dump()
        publish({"type": "response", "role": role, "iteration": iteration,
                 "content": response, "usage": usage, "elapsed": time.perf_counter() - started})
        return response

    try:
        consultant = call("Consultant", settings.consultant_model, prompts["sys_msg"],
                          f"[Problem Statement]:\n{problem}", "Consultant_module", 0)
        plans = "\n\n".join(f"{i+1}. {p['solu_name']}:\n{p['content']}" for i, p in enumerate(consultant["solution_plans"]))
        previous = None
        runtime = None
        review = None
        for iteration in range(settings.revisions + 1):
            if iteration == 0:
                system = prompts["sys_msg_coder"] + prompts["add_first_sys_msg"]
                user = f"[Problem Statement]:\n{problem}\n[Expanded Problem]:\n{consultant['expanded_prob']}\n[Solution Plans]:\n{plans}"
            else:
                system = prompts["sys_msg_coder"] + prompts["add_rev_sys_msg"]
                user = f"[Problem Statement]:\n{problem}\n[Technical Details]:\n{previous['technical_spec']}\n[Previous Code]:\n{previous['python_codes']}\n[Runtime Output]:\n{runtime}\n[Review Comments] need addressed:\n{json.dumps(review)}"
            previous = call("Programmer" if iteration == 0 else "Programmer Revision", settings.programmer_model,
                            system, user, "Progarmmer_module", iteration)
            blocks = python_blocks(previous["python_codes"])
            code = blocks[0] if blocks else previous["python_codes"]
            publish({"type": "start", "role": "Executor", "iteration": iteration})
            try:
                execution = execute_script(code, settings.execution_timeout, generated=True,
                    on_progress=lambda out, err: emit({"type": "progress", "role": "Executor", "iteration": iteration,
                                                       "stdout": out, "stderr": err}))
                runtime = f"Return code: {execution.returncode}\nTimed out: {execution.timed_out}\nSTDOUT:\n{execution.stdout}\nSTDERR:\n{execution.stderr}"
                publish({"type": "execution", "role": "Executor", "iteration": iteration, "result": execution})
            except (ValueError, SyntaxError) as error:
                execution = ExecutionResult("", str(error), 1, 0.0, False)
                runtime = f"Execution validation failed: {error}"
                publish({"type": "execution", "role": "Executor", "iteration": iteration, "result": execution})
            image_review = ""
            if settings.image_review and execution.images:
                import base64
                image_messages = [{"type": "text", "text": f"[Problem statement]:\n{problem}\n[Submitted Code]:\n{code}\n[Runtime Output]:\n{runtime}"}]
                for name, data in list(execution.images.items())[:4]:
                    # GIFs are animations for humans; use only static plot formats with the vision API.
                    if Path(name).suffix.lower() not in {".png", ".jpg", ".jpeg", ".webp"}:
                        continue
                    mime = "image/png" if name.lower().endswith(".png") else "image/webp" if name.lower().endswith(".webp") else "image/jpeg"
                    image_messages.extend([{"type": "text", "text": f"Image filename: {name}"},
                                           {"type": "image_url", "image_url": {"url": f"data:{mime};base64,{base64.b64encode(data).decode()}"}}])
                if len(image_messages) > 1:
                    vision = call("Image Reviewer", settings.reviewer_model, prompts["image_rev_sys_msg"],
                                  image_messages, "Reivew_multimodal_image", iteration)
                    image_review = json.dumps(vision)
            review = call("Reviewer", settings.reviewer_model, prompts["sys_msg_reviewer"],
                          f"[Problem Statement]:\n{problem}\n[Expanded Problem]:\n{consultant['expanded_prob']}\n[Code Submitted]:\n{code}\n[Technical Details]:\n{previous['technical_spec']}\n[Runtime Output]:\n{runtime}\n[Image-based Review]:\n{image_review}",
                          "Reviewer_module", iteration)
            if review["review_decision"] == "accept":
                publish({"type": "complete", "decision": "accept", "rounds": iteration + 1})
                break
        else:
            publish({"type": "complete", "decision": "revision budget reached", "rounds": settings.revisions + 1})
    finally:
        client.close()
    return events
