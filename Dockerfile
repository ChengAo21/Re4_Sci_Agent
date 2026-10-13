FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    MPLBACKEND=Agg \
    RE4_ENABLE_LIVE_EXECUTION=0

RUN useradd --create-home --uid 1000 re4
WORKDIR /app

COPY requirements.txt requirements-live.txt ./
ARG RE4_REQUIREMENTS=requirements.txt
RUN python -m pip install -r "${RE4_REQUIREMENTS}"

COPY --chown=re4:re4 . .
RUN chown re4:re4 /app
USER re4

EXPOSE 7860
HEALTHCHECK --interval=30s --timeout=5s --start-period=30s \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:7860/_stcore/health', timeout=3)"

CMD ["python", "-m", "streamlit", "run", "streamlit_app.py", "--server.address=0.0.0.0", "--server.port=7860", "--server.headless=true", "--browser.gatherUsageStats=false"]
