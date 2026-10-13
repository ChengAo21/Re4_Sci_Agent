"""Re4 Lab — a visual companion to the original scientific-computing agent."""
from __future__ import annotations

import base64
from html import escape
import json
import os
from pathlib import Path
import re

import streamlit as st

from re4_ui.archive import ROOT, archive_choices, code_choices, read_archive, read_trusted_code, redact, solution_plans
from re4_ui.execution import ExecutionResult, execute_script, live_enabled

st.set_page_config(page_title="Re4 Lab · Scientific Computing Agent", page_icon="✳", layout="wide",
                   initial_sidebar_state="collapsed")

STYLE = """
<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&family=Manrope:wght@400;500;600;700;800&display=swap');
html, body, .stApp {font-family: 'DM Sans', sans-serif;}
.stApp {background: #F7F8F5; color: #152D3B;}
.block-container {max-width: 1320px; padding-top: 2.1rem; padding-bottom: 4rem;}
h1,h2,h3 {font-family: 'Manrope','DM Sans',sans-serif; color: #152D3B; letter-spacing:-.035em;}
h1 {font-weight:800!important;} h2 {font-size: 1.8rem!important;} h3 {font-size:1.13rem!important;}
[data-testid="stHeader"] {background: rgba(247,248,245,.92);}
.brand {display:flex; justify-content:space-between; align-items:center; padding-bottom:1.2rem; border-bottom:1px solid #DCE4DF; margin-bottom:1.7rem;}
.brand-name {font-family:Manrope,sans-serif; font-weight:800; font-size:1.55rem; letter-spacing:-.06em;}
.brand-name span {color:#087F8C;} .brand-sub {font-size:.72rem; font-weight:600; letter-spacing:.13em; color:#687F88; text-transform:uppercase;}
.edition {border:1px solid #D7E1DB; padding:.38rem .75rem; border-radius:100px; font-size:.72rem; color:#45626F;}
.eyebrow {font-size:.69rem; letter-spacing:.16em; font-weight:700; text-transform:uppercase; color:#087F8C; margin-bottom:.75rem;}
.hero {padding:1rem 0 1.4rem; max-width:950px;}
.hero h1 {font-size:3.3rem; line-height:1.1; margin:0 0 1rem; max-width:850px;}
.hero p {font-size:1.04rem; line-height:1.65; color:#55707A; max-width:760px;}
.pill {display:inline-block; border-radius:99px; background:#E4F2ED; color:#08766D; font-size:.67rem; font-weight:700; letter-spacing:.06em; padding:.28rem .65rem; margin:0 .3rem .4rem 0;}
.pill-amber {background:#FFF0D8; color:#92610D;} .pill-navy {background:#E8EDF2; color:#3D566B;}
.stage {background:#fff; border:1px solid #DFE7E0; border-radius:16px; padding:1rem 1.1rem; min-height:135px;}
.stage .num {font-size:.68rem; color:#91A59D; font-weight:700; letter-spacing:.08em;}
.stage h3 {margin:.4rem 0 .35rem; font-size:1rem!important; letter-spacing:-.02em;}
.stage p {font-size:.79rem; color:#68808A; line-height:1.5; margin:0;}
.section-intro {color:#607A83; font-size:.92rem; margin-bottom:1.5rem; max-width:830px;}
.panel-title {font-family:Manrope,sans-serif; font-weight:700; font-size:1.25rem; letter-spacing:-.03em; margin-bottom:.6rem;}
.feature-note {font-size:.85rem; color:#55707A; line-height:1.7;}
.meta-line {font-size:.74rem; color:#7C908A; margin-top:.5rem;}
.footer {border-top:1px solid #DFE7E0; margin-top:3rem; padding-top:1rem; color:#789086; font-size:.76rem; display:flex;justify-content:space-between;}
[data-testid="stMetric"] {background:#fff; border:1px solid #DFE7E0; border-radius:13px; padding:.8rem 1rem;}
[data-testid="stMetricValue"] {font-family:Manrope,sans-serif; font-size:1.55rem; font-weight:700;}
[data-testid="stMetricLabel"] {font-size:.75rem; color:#657F88;}
div[data-testid="stVerticalBlockBorderWrapper"] > div {border-radius:15px!important;}
.stButton button, .stDownloadButton button {border-radius:9px; font-weight:600;}
.stCode {font-size:.77rem;} [data-testid="stChatMessage"] {background:#fff; border:1px solid #DFE7E0; border-radius:14px;}
[data-testid="stRadio"] > div {gap:1rem;}
@media(max-width:740px) {.hero h1 {font-size:2.2rem;} .brand-sub {display:none;} .edition{font-size:.62rem;} .block-container{padding-top:1.2rem;}}
</style>
"""
st.markdown(STYLE, unsafe_allow_html=True)


def local_file(relative: str) -> Path | None:
    path = (ROOT / relative).resolve()
    if path.is_relative_to(ROOT) and path.is_file():
        return path
    return None


@st.cache_data(show_spinner=False, ttl=30)
def load_manifest() -> dict:
    path = ROOT / "assets" / "results" / "manifest.json"
    if path.is_file():
        try:
            return json.loads(path.read_text(encoding="utf-8-sig"))
        except (OSError, ValueError):
            return {"results": []}
    return {"results": []}


def show_image(path: Path | bytes, *, caption: str | None = None) -> None:
    if isinstance(path, Path) and path.suffix.lower() == ".gif":
        encoded = base64.b64encode(path.read_bytes()).decode()
        st.markdown(f'<img src="data:image/gif;base64,{encoded}" alt="Cylinder flow animation" style="width:100%;border-radius:12px;">', unsafe_allow_html=True)
        if caption:
            st.caption(caption)
    else:
        st.image(str(path) if isinstance(path, Path) else path, caption=caption, width="stretch")


def stage_strip() -> None:
    stages = [("01", "Rewriting", "Consultant", "Expand the problem. Compare solution plans."),
              ("02", "Resolution", "Programmer + Python", "Choose a method. Write and execute a solver."),
              ("03", "Review", "Reviewer", "Inspect the code, logs and numerical results."),
              ("04", "Revision", "Programmer + Python", "Use concrete feedback to improve the next attempt.")]
    for col, (num, title, role, detail) in zip(st.columns(4), stages):
        with col:
            st.markdown(f'<div class="stage"><div class="num">{num} · {role}</div><h3>{title}</h3><p>{detail}</p></div>', unsafe_allow_html=True)


@st.cache_data(show_spinner=False)
def cylinder_states(path: str) -> dict:
    import numpy as np
    with np.load(path, allow_pickle=False) as data:
        return {name: data[name] for name in ("times", "u", "v", "x", "y")}


def explore_cylinder() -> None:
    path = local_file("assets/results/cylinder_snapshots.npz")
    if not path:
        return
    # A checkbox avoids loading a 12 MB scientific field cache on every gallery visit.
    with st.expander("Explore the actual saved time steps"):
        enabled = st.checkbox("Load interactive cylinder fields", key="load_cylinder_fields")
        st.caption("The slider selects saved solver states. No interpolated or synthetic frames are added.")
        if not enabled:
            return
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle
        data = cylinder_states(str(path))
        times = [float(t) for t in data["times"]]
        time_value = st.select_slider("Simulation time", options=times, value=times[-1], format_func=lambda v: f"t = {v:.3f}")
        index = times.index(time_value)
        speed = np.hypot(data["u"][index], data["v"][index])
        figure, axis = plt.subplots(figsize=(11, 4.2))
        figure.set_facecolor("#F7F8F5")
        plot = axis.pcolormesh(data["x"], data["y"], speed.T, shading="auto", cmap="RdBu_r", vmin=0, vmax=1.71, rasterized=True)
        axis.add_patch(Circle((0, 0), .5, facecolor="#152D3B", edgecolor="white", linewidth=.7))
        axis.set(xlabel="x / D", ylabel="y / D", title=f"Velocity magnitude · t = {time_value:.3f}")
        axis.set_aspect("equal")
        figure.colorbar(plot, ax=axis, label="|u| / U∞", shrink=.8, pad=.025)
        figure.tight_layout()
        st.pyplot(figure, width="stretch")
        plt.close(figure)
        a, b, c = st.columns(3)
        a.metric("Saved time", f"{time_value:.3f}")
        b.metric("Snapshot", f"{index+1} / {len(times)}")
        c.metric("Grid", f"{len(data['x'])} × {len(data['y'])}")


def result_metadata(result: dict) -> None:
    if result.get("source"):
        st.caption("Archived source")
        st.code(result["source"], language="text", wrap_lines=True)
    if result.get("metrics"):
        st.json(result["metrics"], expanded=True)
    if result.get("description"):
        st.write(result["description"])


def gallery() -> None:
    st.markdown('<div class="hero"><div class="eyebrow">The Re4 research companion</div><h1>From a problem<br>to a working solver.</h1><p>Explore numerical results, follow the agent’s trial-and-error process, and run the archived solvers. Four stages. Three collaborating roles. Every result has a source.</p></div>', unsafe_allow_html=True)
    manifest = load_manifest()
    results = manifest.get("results", [])
    stats = st.columns(3)
    stats[0].metric("Operational stages", "4")
    stats[1].metric("Collaborating LLM roles", "3")
    stats[2].metric("Displayed numerical cases", len(results))
    st.write("")
    stage_strip()
    st.write("")
    cylinder = next((r for r in results if "cylinder" in r.get("id", "").lower()), None)
    if cylinder:
        st.markdown("## A solver in motion")
        with st.container(border=True):
            image_col, text_col = st.columns([1.8, 1], gap="large")
            with image_col:
                image = local_file(cylinder.get("gif") or cylinder.get("image", ""))
                if image:
                    show_image(image)
            with text_col:
                st.markdown('<span class="pill">FEATURED CASE</span><span class="pill pill-amber">ARCHIVED CODE</span>', unsafe_allow_html=True)
                st.markdown(f"### {cylinder.get('title', 'Flow past a cylinder')}")
                st.markdown(f'<p class="feature-note">{escape(cylinder.get("subtitle", "Incompressible flow, numerical boundary treatment, and an evolving wake."))}</p>', unsafe_allow_html=True)
                if cylinder.get("description"):
                    st.write(cylinder["description"])
                st.caption("Animation provenance and run settings are recorded with the result below.")
                with st.expander("Run settings & source"):
                    result_metadata(cylinder)
        explore_cylinder()
    st.markdown("## Numerical gallery")
    st.markdown('<p class="section-intro">Results produced from the bundled solver archives. These are individual demonstrations, not new aggregate benchmark scores.</p>', unsafe_allow_html=True)
    if not results:
        st.info("Gallery assets have not been generated yet. The archived interactions are available in Agent replay.")
        return
    categories = sorted({r.get("category", "Scientific computing") for r in results})
    category = st.selectbox("Explore a field", ["All cases"] + categories)
    visible = [r for r in results if category == "All cases" or r.get("category", "Scientific computing") == category]
    for start in range(0, len(visible), 2):
        for col, result in zip(st.columns(2, gap="large"), visible[start:start+2]):
            with col:
                with st.container(border=True):
                    st.markdown(f'<span class="pill pill-navy">{escape(result.get("category", "NUMERICAL RESULT"))}</span>', unsafe_allow_html=True)
                    st.markdown(f"### {result.get('title', result.get('id', 'Result'))}")
                    if result.get("subtitle"):
                        st.caption(result["subtitle"])
                    image = local_file(result.get("image", ""))
                    if image:
                        show_image(image)
                    else:
                        st.info("No static plot is available for this case.")
                    with st.expander("Evidence, settings & additional plots"):
                        result_metadata(result)
                        for item in result.get("images", []):
                            relative = item if isinstance(item, str) else item.get("path", "")
                            other = local_file(relative)
                            if other and other != image:
                                show_image(other, caption=None if isinstance(item, str) else item.get("caption"))
                        if image:
                            st.download_button("Download figure", image.read_bytes(), image.name,
                                               key=f"download_{result.get('id', start)}")
    if manifest.get("generated_at"):
        st.caption(f"Gallery generated: {manifest['generated_at']} · Provenance: docs/results_provenance.md")


def show_execution(result: ExecutionResult, key: str) -> None:
    a, b, c = st.columns(3)
    a.metric("Wall time", f"{result.elapsed:.2f} s")
    b.metric("Process exit", "Timeout" if result.timed_out else str(result.returncode))
    c.metric("Generated plots", len(result.images))
    if result.timed_out:
        st.warning("The process reached the configured wall-clock limit and was stopped. Partial outputs are retained.")
    elif result.succeeded:
        st.success("Python execution completed. A zero exit code alone does not establish numerical accuracy.")
    else:
        st.error("Python execution failed. Inspect the actual traceback below.")
    log_tab, figure_tab = st.tabs(["Console", "Generated figures"])
    with log_tab:
        st.code(result.stdout or "No standard output.", language="text", wrap_lines=True)
        if result.stderr:
            with st.expander("Standard error / warnings", expanded=not result.succeeded):
                st.code(result.stderr, language="text", wrap_lines=True)
    with figure_tab:
        if not result.images:
            st.info("This execution did not save any supported images.")
        for name, data in result.images.items():
            show_image(data, caption=name)
    if result.bundle:
        st.download_button("Download this run · code, logs and outputs", result.bundle,
                           "re4_solver_run.zip", "application/zip", key=f"run_download_{key}")


def replay() -> None:
    st.markdown("## Inside the agent loop")
    st.markdown('<p class="section-intro">Read the original Consultant, Programmer and Reviewer exchanges, one attempt at a time. This is a replay of saved reports; browsing makes no API calls.</p>', unsafe_allow_html=True)
    choices = archive_choices()
    if not choices:
        st.warning("No archived Markdown reports were found.")
        return
    labels = list(choices)
    default = next((i for i, name in enumerate(labels) if "Cylinder flow" in name and "3_" in name), 0)
    selection = st.selectbox("Archived experiment", labels, index=default)
    archive = read_archive(choices[selection])
    st.markdown('<span class="pill pill-amber">ARCHIVED INTERACTION</span><span class="pill pill-navy">NO API CALLS</span>', unsafe_allow_html=True)
    info = st.columns(3)
    info[0].metric("Saved code attempts", len(archive.rounds))
    info[1].metric("Saved final review", archive.final_decision.capitalize())
    info[2].metric("Source format", "Markdown report")
    with st.expander("Original task & report provenance"):
        st.markdown(archive.problem)
        st.code(archive.path.relative_to(ROOT).as_posix(), language="text", wrap_lines=True)
        st.download_button("Download archived interaction", archive.raw, archive.path.name, "text/markdown")
    if not archive.rounds:
        st.warning("No Python code attempts were found in this archived report.")
        return
    consultant_tab, attempts_tab = st.tabs(["01 · Consultant / Rewriting", "02–04 · Resolve, review, revise"])
    with consultant_tab:
        with st.chat_message("assistant", avatar="🧠"):
            st.markdown("**Consultant · expanded problem**")
            st.markdown(archive.expanded_problem or "Not recorded in this report.")
        with st.chat_message("assistant", avatar="🧭"):
            st.markdown("**Consultant · alternative solution plans**")
            plans = solution_plans(archive.plans)
            if plans:
                for i, (title, content) in enumerate(plans):
                    with st.expander(f"Plan {i+1} · {title}", expanded=i == 0):
                        st.markdown(content)
            else:
                st.markdown(archive.plans or "Not recorded in this report; the Consultant may have been disabled for this trial.")
    with attempts_tab:
        labels = ["Initial attempt" if i == 0 else f"Revision {i}" for i in range(len(archive.rounds))]
        selected = st.radio("Follow the trajectory", range(len(labels)), format_func=lambda i: labels[i], horizontal=True)
        attempt = archive.rounds[selected]
        st.caption(f"{labels[selected]} → Python execution → Reviewer feedback")
        with st.chat_message("assistant", avatar="💻"):
            st.markdown(f"**Programmer · {labels[selected].lower()}**")
            if attempt.specification:
                st.markdown(attempt.specification)
            with st.expander("Inspect the generated Python", expanded=False):
                st.code(attempt.code, language="python", line_numbers=True)
                st.download_button("Download this script", attempt.code, f"re4_attempt_{selected}.py", "text/x-python", key=f"code_{selection}_{selected}")
        with st.chat_message("assistant", avatar="⚙️"):
            st.markdown("**Executor · archived runtime output**")
            st.code(attempt.runtime or "No runtime output was saved for this attempt.", language="text", wrap_lines=True)
        with st.chat_message("assistant", avatar="🔎"):
            st.markdown("**Reviewer · archived diagnostic feedback**")
            st.markdown(attempt.review or "No review comments were saved for this attempt.")
            if attempt.image_review:
                with st.expander("Archived image-review notes"):
                    st.markdown(attempt.image_review)
        st.markdown("### Re-run this archived attempt")
        st.caption("Runs the saved script in a fresh workspace with a time limit. Open Matplotlib figures are automatically exported. The historical review above is not a review of this new execution.")
        timeout = st.select_slider("Execution limit", options=[15, 30, 60, 120, 180], value=60, key="replay_timeout")
        run_key = f"{archive.path.as_posix()}:{selected}"
        if st.button("Run archived Python", type="primary", key="replay_run"):
            try:
                with st.status("Executing archived Python…", expanded=True) as status:
                    st.write("Created a clean workspace and copied the bundled data files.")
                    console = st.empty()
                    def update_console(out: str, err: str) -> None:
                        console.code(out or err or "Python process started…", language="text", wrap_lines=True)
                    st.session_state["replay_result"] = (run_key, execute_script(attempt.code, timeout, on_progress=update_console))
                    status.update(label="Execution finished", state="complete")
            except (ValueError, SyntaxError, OSError) as error:
                st.error(f"Execution could not start: {redact(str(error))}")
        saved = st.session_state.get("replay_result")
        if saved and saved[0] == run_key:
            show_execution(saved[1], "replay")


def render_live_content(role: str, content: dict) -> None:
    if role == "Consultant":
        st.markdown(content.get("expanded_prob", ""))
        for i, plan in enumerate(content.get("solution_plans", [])):
            with st.expander(f"Plan {i+1} · {plan.get('solu_name', '')}"):
                st.markdown(plan.get("content", ""))
    elif role.startswith("Programmer"):
        st.markdown(content.get("technical_spec", ""))
        with st.expander("Generated Python"):
            from re4_ui.archive import python_blocks
            blocks = python_blocks(content.get("python_codes", ""))
            st.code(blocks[0] if blocks else content.get("python_codes", ""), language="python")
    elif role == "Reviewer":
        decision = content.get("review_decision", "")
        st.markdown(f"**Decision: {decision.upper()}**")
        st.markdown(content.get("review_summary", ""))
        for comment in content.get("review_comments", []):
            with st.expander(f"{comment.get('severity', '')} · {comment.get('category', '')} · {comment.get('issue', '')}"):
                st.markdown(comment.get("recommendation", ""))
    else:
        st.json(content)


def live_lab() -> None:
    st.markdown("## Run lab")
    st.markdown('<p class="section-intro">Reproduce a bundled final-stage solver, or connect your own API to watch Re4 generate, execute, review and revise code.</p>', unsafe_allow_html=True)
    trusted_tab, live_tab = st.tabs(["Archived solvers · no API", "Live agent · your API"])
    with trusted_tab:
        choices = code_choices()
        if choices:
            labels = list(choices)
            default = next((i for i, name in enumerate(labels) if "SodShock" in name and "GPT_ans" in name and "sample1" in name), 0)
            choice = st.selectbox("Final-stage archived code", labels, index=default)
            path, index = choices[choice]
            code = read_trusted_code(path, index)
            with st.expander("Source Python & provenance"):
                st.code(path.relative_to(ROOT).as_posix(), language="text", wrap_lines=True)
                st.code(code, language="python", line_numbers=True)
            timeout = st.select_slider("Wall-clock limit", options=[15, 30, 60, 120, 180], value=60)
            if st.button("Execute selected solver", type="primary"):
                try:
                    with st.status("Running the selected archived solver…", expanded=True) as status:
                        st.write("The solver will save figures and logs in a temporary workspace. Open Matplotlib figures are automatically exported.")
                        console = st.empty()
                        def update_console(out: str, err: str) -> None:
                            console.code(out or err or "Python process started…", language="text", wrap_lines=True)
                        st.session_state["trusted_result"] = (choice, execute_script(code, timeout, on_progress=update_console))
                        status.update(label="Archived solver finished", state="complete")
                except (ValueError, SyntaxError, OSError) as error:
                    st.error(f"Execution could not start: {redact(str(error))}")
            saved = st.session_state.get("trusted_result")
            if saved and saved[0] == choice:
                show_execution(saved[1], "trusted")
        else:
            st.info("No archived final-stage code was found.")
    with live_tab:
        if not live_enabled():
            st.info("Live agent execution is disabled on this host. The gallery, interaction replay and archived solver runs are available without an API key.")
            st.caption("On a controlled local host, enable RE4_ENABLE_LIVE_EXECUTION=1 before launching Streamlit to activate this panel.")
            return
        from re4_ui.live import LiveSettings, problems, run_pipeline
        presets = problems()
        names = list(presets)
        task_name = st.selectbox("Original task prompt", names, index=names.index("burgers_problem") if "burgers_problem" in names else 0)
        with st.form("live_settings"):
            problem = st.text_area("Problem statement", presets[task_name], height=230)
            st.markdown("#### API connection")
            key = st.text_input("API key", type="password", help="Used in this session; never included in exported logs.")
            base_url = st.text_input("OpenAI-compatible API base URL", value=os.getenv("RE4_API_BASE", "https://api.openai.com/v1"))
            model_default = os.getenv("RE4_MODEL", "")
            c1, c2, c3 = st.columns(3)
            cons = c1.text_input("Consultant model", value=os.getenv("RE4_CONSULTANT_MODEL", model_default))
            prog = c2.text_input("Programmer model", value=os.getenv("RE4_PROGRAMMER_MODEL", model_default))
            rev = c3.text_input("Reviewer model", value=os.getenv("RE4_REVIEWER_MODEL", model_default))
            a, b, c = st.columns(3)
            revisions = a.selectbox("Maximum revisions", [0, 1, 2], index=2)
            timeout = b.selectbox("Python timeout (seconds)",[30, 60, 120, 180, 300, 600, 1200, 1800],index=5,)
            token_limit = c.selectbox("Output-token limit per call", [4000, 8000, 12000, 16000], index=1)
            vision = st.checkbox("Enable image review (requires a vision-capable Reviewer model)")
            approved = st.checkbox("I authorize calls to my API provider and execution of the generated Python on this controlled host. Provider charges may apply.")
            submit = st.form_submit_button("Start a real Re4 run", type="primary")
        if submit:
            api_key = key.strip() or os.getenv("OPENAI_API_KEY", "")
            if not approved or not api_key or not all([cons.strip(), prog.strip(), rev.strip(), problem.strip()]):
                st.error("Provide a task, API key and all role models, then authorize the run.")
                return
            if not base_url.startswith("https://") and not base_url.startswith(("http://localhost", "http://127.0.0.1")):
                st.error("Use an HTTPS API endpoint, or a localhost endpoint for a local model server.")
                return
            settings = LiveSettings(api_key, base_url, cons.strip(), prog.strip(), rev.strip(), revisions, timeout, token_limit, vision)
            st.session_state["live_events"] = []
            st.markdown("### Live trajectory")
            event_area = st.container()
            current = {}
            status = st.status("Re4 is running…", expanded=True)
            def emit(event: dict) -> None:
                role = event.get("role", "")
                iteration = event.get("iteration", 0)
                with event_area:
                    if event["type"] == "start":
                        status.write(f"{role} · attempt {iteration + 1}" + (f" · {event['model']}" if event.get("model") else ""))
                        with st.chat_message("assistant", avatar="🧠" if role == "Consultant" else "🔎" if "Reviewer" in role else "💻"):
                            st.markdown(f"**{role} · attempt {iteration+1}**")
                            current[(role, iteration)] = st.empty()
                    elif event["type"] == "token":
                        current[(role, iteration)].code(event["text"], language="json", wrap_lines=True)
                    elif event["type"] == "response":
                        placeholder = current[(role, iteration)]
                        with placeholder.container():
                            render_live_content(role, event["content"])
                            total = (event.get("usage") or {}).get("total_tokens")
                            st.caption(f"{event['elapsed']:.1f} s" + (f" · {total:,} tokens" if total else " · token usage not reported by provider"))
                    elif event["type"] == "execution":
                        with current[(role, iteration)].container():
                            show_execution(event["result"], f"live_{iteration}")
                    elif event["type"] == "progress":
                        current[(role, iteration)].code(event["stdout"] or event["stderr"] or "Python process started…", language="text", wrap_lines=True)
                    elif event["type"] == "complete":
                        status.update(label=f"Run finished · {event['decision']} · {event['rounds']} attempts", state="complete")
                if event["type"] not in {"token", "progress"}:
                    st.session_state["live_events"].append(event)
            try:
                run_pipeline(problem, settings, emit)
            except Exception as error:
                status.update(label="Run stopped", state="error")
                message = redact(str(error)).replace(api_key, "[REDACTED]")
                st.error(f"The provider or pipeline returned an error: {message[:1500]}")
            # The key remains an input widget value, never an event or downloadable setting.
        elif st.session_state.get("live_events"):
            st.markdown("### Previous live run in this session")
            for i, event in enumerate(st.session_state["live_events"]):
                if event["type"] == "response":
                    with st.expander(f"{event['role']} · attempt {event['iteration'] + 1}"):
                        render_live_content(event["role"], event["content"])
                elif event["type"] == "execution":
                    with st.expander(f"Executor · attempt {event['iteration'] + 1}"):
                        show_execution(event["result"], f"previous_{i}")


def about() -> None:
    st.markdown("## Scientific computing, with a review loop")
    st.markdown("Re4 connects three language-model roles with a Python execution tool. Its four operational stages are **Rewriting, Resolution, Review and Revision**.")
    stage_strip()
    st.write("")
    left, right = st.columns([1.3, 1], gap="large")
    with left:
        if (ROOT / "schematic_re4.png").is_file():
            st.image(str(ROOT / "schematic_re4.png"), width="stretch")
    with right:
        st.markdown("### What this companion shows")
        st.markdown("- **Results gallery:** numerical plots and animations with source records.\n- **Agent replay:** original text, generated Python, runtime reports and review feedback.\n- **Run lab:** fresh execution of archived scripts, plus an opt-in connection to your LLM provider.")
        st.markdown("### Reading the evidence")
        st.write("A script finishing successfully is different from a numerically validated solution. The app preserves failures and historical reviewer comments so that implementation choices can be examined directly.")
        st.caption("The original code remains in src/pycodes_ver and src/ipynb_ver. The presentation adapter adds a UI without making API calls on import.")


st.markdown('<div class="brand"><div><div class="brand-name">Re<span>4</span> <span style="color:#607A83;font-weight:500">Lab</span></div><div class="brand-sub">Scientific computing agent</div></div><div class="edition">REWRITING · RESOLUTION · REVIEW · REVISION</div></div>', unsafe_allow_html=True)
page = st.radio("Navigation", ["Results gallery", "Agent replay", "Run lab", "About Re4"], horizontal=True, label_visibility="collapsed")
st.write("")
if page == "Results gallery":
    gallery()
elif page == "Agent replay":
    replay()
elif page == "Run lab":
    live_lab()
else:
    about()
st.markdown('<div class="footer"><span>Re4 · Rewriting, Resolution, Review and Revision</span><span>Research code · Traceable results</span></div>', unsafe_allow_html=True)
