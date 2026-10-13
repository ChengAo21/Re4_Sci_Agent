from pathlib import Path

from streamlit.testing.v1 import AppTest

from re4_ui.archive import ROOT


def test_all_pages_load_without_errors(monkeypatch):
    monkeypatch.delenv("RE4_ENABLE_LIVE_EXECUTION", raising=False)
    app = AppTest.from_file(str(ROOT / "streamlit_app.py")).run(timeout=30)
    assert not app.exception
    assert any(m.label == "Operational stages" and m.value == "4" for m in app.metric)
    for page in ["Agent replay", "Run lab", "About Re4"]:
        app.radio[0].set_value(page).run(timeout=30)
        assert not app.exception, page


def test_run_lab_executes_actual_archived_sod(monkeypatch):
    monkeypatch.delenv("RE4_ENABLE_LIVE_EXECUTION", raising=False)
    app = AppTest.from_file(str(ROOT / "streamlit_app.py")).run(timeout=30)
    app.radio[0].set_value("Run lab").run(timeout=30)
    app.button[0].click().run(timeout=60)
    assert not app.exception
    assert any(m.label == "Process exit" and m.value == "0" for m in app.metric)
    assert any(m.label == "Generated plots" and int(m.value) > 0 for m in app.metric)
    assert any("Solution summary" in block.value for block in app.code)


def test_true_cylinder_time_slider_loads(monkeypatch):
    monkeypatch.delenv("RE4_ENABLE_LIVE_EXECUTION", raising=False)
    app = AppTest.from_file(str(ROOT / "streamlit_app.py")).run(timeout=30)
    app.checkbox(key="load_cylinder_fields").check().run(timeout=30)
    assert not app.exception
    assert any(m.label == "Saved time" and m.value == "20.000" for m in app.metric)
    slider = next(s for s in app.select_slider if s.label == "Simulation time")
    slider.set_value(0.0).run(timeout=30)
    assert not app.exception
    assert any(m.label == "Saved time" and m.value == "0.000" for m in app.metric)


def test_live_form_is_explicit_and_no_requests_on_load(monkeypatch):
    monkeypatch.setenv("RE4_ENABLE_LIVE_EXECUTION", "1")
    app = AppTest.from_file(str(ROOT / "streamlit_app.py")).run(timeout=30)
    app.radio[0].set_value("Run lab").run(timeout=30)
    assert not app.exception
    assert any(field.label == "API key" and field.value == "" for field in app.text_input)
    assert any("Provider charges may apply" in field.label and not field.value for field in app.checkbox)
