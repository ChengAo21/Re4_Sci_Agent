# Verification — 2026-10-13

- Python 3.12.14, Streamlit 1.65.0.
- 20 automated tests passed in 15.60 seconds, including actual archived execution through Streamlit AppTest.
- All 10 showcase solvers were rerun. Their source hashes and generated images/logs match the manifest.
- Cylinder animation: 66 actual states, original 220×120 grid, t=0–20.
- Browser verification: Gallery, archived interactions, live console, and a fresh Sod execution (exit 0, 2.16 seconds, one generated figure).
- Local health endpoint returned HTTP 200 and `ok`.
- Browser time slider changed between the actual initial and final cylinder states (1/66 and 66/66).
- No broken local documentation links.
- No paid LLM calls were made. The optional live-provider path requires the owner's configured API and was not tested against an external model.
- Docker deployment configuration is included; an image build was not performed because Docker is unavailable on this host.
