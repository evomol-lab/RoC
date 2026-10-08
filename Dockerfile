# The Room of Conformations: web interface + command line, with ClustENM (OpenMM).
#
# Image names are fully qualified so the same file builds with Docker and Podman.
#
#   podman build -t docker.io/evomol/roc:2.0.0 .
#   podman run --rm -p 5050:5050 docker.io/evomol/roc:2.0.0                    # web interface
#   podman run --rm -v "$PWD":/data -w /data docker.io/evomol/roc:2.0.0 \
#       python /opt/roc/roc.py structure.pdb --outdir results                # command line

ARG PYTHON_IMAGE=docker.io/library/python:3.12-slim

FROM ${PYTHON_IMAGE} AS build
ENV PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1
RUN python -m venv /opt/venv
COPY requirements.txt /tmp/requirements.txt
RUN /opt/venv/bin/pip install -r /tmp/requirements.txt openmm==8.6.1 pdbfixer==1.12.0

FROM ${PYTHON_IMAGE}
LABEL org.opencontainers.image.title="The Room of Conformations" \
      org.opencontainers.image.description="All-atom conformational ensembles from elastic network normal modes" \
      org.opencontainers.image.source="https://github.com/jpmslima/RoC" \
      org.opencontainers.image.vendor="EvoMol-Lab, BioME, UFRN" \
      org.opencontainers.image.licenses="MIT"
ENV PATH=/opt/venv/bin:$PATH \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    ROC_HOST=0.0.0.0 \
    ROC_PORT=5050 \
    ROC_JOBS_DIR=/opt/roc/roc_jobs
COPY --from=build /opt/venv /opt/venv
WORKDIR /opt/roc
COPY roc.py RoCGUI.py DOCUMENTATION.md README.md LICENSE.txt ./
COPY templates/ templates/
COPY static/ static/
# Jobs of the web interface; mount a host directory here to keep them.
RUN mkdir -p roc_jobs && chmod 777 roc_jobs
EXPOSE 5050
STOPSIGNAL SIGINT
CMD ["python", "/opt/roc/RoCGUI.py"]
