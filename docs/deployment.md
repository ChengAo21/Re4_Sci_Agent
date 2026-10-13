# Deploy the Streamlit workspace

The app is ready for local execution and Docker-based hosting. The default deployment shows regenerated results and archived interactions without API credentials. This copied working directory has no Git remote configured; the launch and build commands do not initialize or publish a repository.

## Streamlit Community Cloud

1. Upload the prepared repository files to your GitHub repository. Keep `streamlit_app.py`, `requirements.txt`, `re4_ui/`, `assets/`, and `src/` together at the repository root.
2. Sign in to Streamlit Community Cloud and create an app. Select the repository, branch and `streamlit_app.py` entry point; choose Python 3.12 in the deployment settings.
3. Deploy. Gallery, Replay and archived-program execution use `requirements.txt`. No model key is needed for these views.

See the official [deployment guide](https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app). Optional live dependencies are in `requirements-live.txt`; install those only on a private deployment where you intend to execute newly generated code. Configure credentials through the host's secrets interface, never in tracked source files.

## Docker on a server

Build and start from the repository root:

```bash
docker build -t re4-streamlit .
docker run --rm -p 7860:7860 re4-streamlit
```

Open `http://localhost:7860`. The image runs as a non-root user, serves Streamlit on port 7860, and includes a health check at `/_stcore/health`. Runtime work is temporary unless you add persistent storage. The container is an application package; it is not a sandbox for arbitrary generated programs.

For a private live-agent deployment:

```bash
docker build --build-arg RE4_REQUIREMENTS=requirements-live.txt -t re4-streamlit-live .
docker run --rm -p 7860:7860 --env-file .env re4-streamlit-live
```

Set `RE4_ENABLE_LIVE_EXECUTION=1` and the provider configuration in your private `.env` file before running. `.dockerignore` excludes that file from the image. Consult Streamlit's [Docker deployment guide](https://docs.streamlit.io/deploy/tutorials/docker) for host and networking details.

## Hugging Face Docker Space

Create a Space with the **Docker** SDK, then upload this repository. On the Space copy of `README.md`, add the following metadata at the very top (before the GitHub presentation content):

```yaml
---
title: Re4 Scientific Computing
emoji: 🧪
colorFrom: blue
colorTo: green
sdk: docker
app_port: 7860
pinned: false
---
```

The included Dockerfile already listens on port 7860. The default public demonstration requires no secrets. Account access is required to create and deploy the Space; an actual hosted URL is available only after its build and health checks succeed. See the official [Docker Spaces documentation](https://huggingface.co/docs/hub/spaces-sdks-docker).

## Verify a deployment

- The Gallery loads its images and animation without broken paths.
- Replay opens an archived transcript without contacting a model provider.
- A short archived run returns console output and a figure.
- `/_stcore/health` returns `ok`.
- Live execution remains disabled on public demonstrations.

The repository contains no deployment credentials and no assertion that a public instance is already running.
