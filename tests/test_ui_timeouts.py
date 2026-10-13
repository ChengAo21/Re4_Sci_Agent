"""Exercise timeout choices through the rendered UI without API calls or long runs."""
from __future__ import annotations

import pytest
from streamlit.testing.v1 import AppTest

from re4_ui.archive import ROOT
from re4_ui import execution, live
from re4_ui.execution import (
    DEFAULT_EXECUTION_TIMEOUT,
    EXECUTION_TIMEOUT_OPTIONS,
    MAX_EXECUTION_TIMEOUT,
    ExecutionResult,
)


def open_page(monkeypatch, page: str):
    monkeypatch.setenv("RE4_ENABLE_LIVE_EXECUTION", "1")
    app = AppTest.from_file(str(ROOT / "streamlit_app.py")).run(timeout=30)
    assert not app.exception
    app.radio[0].set_value(page).run(timeout=30)
    assert not app.exception
    return app


def button_with_label(app, label: str):
    return next(button for button in app.button if button.label == label)


@pytest.mark.parametrize(
    "page, kind, key",
    [
        ("Agent replay", "select_slider", "replay_timeout"),
        ("Run lab", "select_slider", "trusted_timeout"),
        ("Run lab", "selectbox", "live_timeout"),
    ],
)
def test_all_execution_views_offer_the_same_ten_minute_default_and_thirty_minute_limit(
    monkeypatch, page, kind, key
):
    app = open_page(monkeypatch, page)
    control = getattr(app, kind)(key=key)
    assert DEFAULT_EXECUTION_TIMEOUT == 600
    assert MAX_EXECUTION_TIMEOUT == 1800
    assert tuple(int(option) for option in control.options) == tuple(EXECUTION_TIMEOUT_OPTIONS)
    assert control.value == DEFAULT_EXECUTION_TIMEOUT

    control.set_value(MAX_EXECUTION_TIMEOUT)
    if key == "live_timeout":
        # Streamlit forms send their widget values only on submission. Missing
        # credentials/authorization keep this submission from starting a run.
        button_with_label(app, "Start a real Re4 run").click().run(timeout=30)
    else:
        app.run(timeout=30)
    assert not app.exception
    assert getattr(app, kind)(key=key).value == MAX_EXECUTION_TIMEOUT


@pytest.mark.parametrize(
    "page, key, button_label",
    [
        ("Agent replay", "replay_timeout", "Run archived Python"),
        ("Run lab", "trusted_timeout", "Execute selected solver"),
    ],
)
def test_archived_execution_receives_selected_limit_and_hides_results_for_another_limit(
    monkeypatch, page, key, button_label
):
    calls = []
    marker = "TIMEOUT_TEST_EXECUTED_WITH_1800_SECONDS"

    def fake_execute(code, timeout, **kwargs):
        calls.append((code, timeout))
        return ExecutionResult(marker, "", 0, 0.01, False)

    monkeypatch.setattr(execution, "execute_script", fake_execute)
    app = open_page(monkeypatch, page)
    app.select_slider(key=key).set_value(1800).run(timeout=30)
    button_with_label(app, button_label).click().run(timeout=30)
    assert not app.exception
    assert len(calls) == 1
    assert calls[0][0].strip()
    assert calls[0][1] == 1800
    assert any(block.value == marker for block in app.code)
    assert any(metric.label == "Process exit" and metric.value == "0" for metric in app.metric)

    # Selecting another time budget must not present the previous run as its result.
    app.select_slider(key=key).set_value(600).run(timeout=30)
    assert not app.exception
    assert len(calls) == 1  # Changing a setting does not execute a script automatically.
    assert not any(block.value == marker for block in app.code)
    assert not any(metric.label == "Process exit" for metric in app.metric)


def test_live_form_passes_selected_thirty_minute_limit_to_pipeline_without_api_calls(monkeypatch):
    calls = []

    def fake_pipeline(problem, settings, emit):
        calls.append((problem, settings))
        emit({"type": "complete", "decision": "accept", "rounds": 1})
        return []

    monkeypatch.setattr(live, "run_pipeline", fake_pipeline)
    app = open_page(monkeypatch, "Run lab")
    app.selectbox(key="live_timeout").set_value(1800)
    values = {
        "API key": "test-key-with-no-provider-access",
        "OpenAI-compatible API base URL": "https://example.invalid/v1",
        "Consultant model": "test-model",
        "Programmer model": "test-model",
        "Reviewer model": "test-model",
    }
    for field in app.text_input:
        if field.label in values:
            field.set_value(values[field.label])
    authorization = next(field for field in app.checkbox if "Provider charges may apply" in field.label)
    authorization.check()
    button_with_label(app, "Start a real Re4 run").click().run(timeout=30)

    assert not app.exception
    assert len(calls) == 1
    assert calls[0][0].strip()
    assert calls[0][1].execution_timeout == 1800
    assert calls[0][1].api_key == values["API key"]


@pytest.mark.parametrize("requested_timeout, effective_timeout", [(600, 600), (3600, 1800)])
def test_backend_preserves_a_budget_above_three_minutes_without_waiting(
    monkeypatch, requested_timeout, effective_timeout
):
    clock = {"time": 0.0}

    class SimulatedProcess:
        returncode = None

        def __init__(self, *args, **kwargs):
            launcher = (kwargs["cwd"] / "_launch.py").read_text(encoding="utf-8")
            assert f"resource.setrlimit(resource.RLIMIT_CPU, ({effective_timeout}, {effective_timeout + 1}))" in launcher

        def poll(self):
            if self.returncode is None and clock["time"] >= 400:
                self.returncode = 0
            return self.returncode

        def kill(self):
            self.returncode = -9

        def wait(self, timeout):
            return self.returncode

    monkeypatch.setattr(execution.subprocess, "Popen", SimulatedProcess)
    monkeypatch.setattr(execution.time, "perf_counter", lambda: clock["time"])
    monkeypatch.setattr(execution.time, "sleep", lambda _: clock.update(time=clock["time"] + 200))

    # The process remains running at 200 s and finishes at 400 s. A hidden 180 s
    # clamp would kill it at the second poll. Neither a child nor a sleep occurs.
    result = execution.execute_script("print('simulated numerical solver')", timeout=requested_timeout)
    assert result.succeeded
    assert not result.timed_out
    assert result.elapsed == 400
