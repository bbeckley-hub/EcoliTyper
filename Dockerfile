# syntax=docker/dockerfile:1.6

FROM mambaorg/micromamba:1.5.8-jammy

LABEL maintainer="Brown Beckley <brownbeckley94@gmail.com>"
LABEL description="EcoliTyper - Complete E. coli Typing Pipeline"

USER root

# --- System dependencies ------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
        jq \
        git \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/ecolityper

# --- Copy environment spec FIRST so code edits don't invalidate the env layer -
COPY environment.yml /opt/ecolityper/environment.yml

# --- Create the env with a BuildKit cache mount ------------------------------
# No `micromamba clean` here — /opt/conda/pkgs IS the cache and must be kept
# for retry-safe rebuilds. The cache mount survives --no-cache.
RUN --mount=type=cache,target=/opt/conda/pkgs \
    micromamba create -y -n ecolityper -f /opt/ecolityper/environment.yml

# --- Copy the project ---------------------------------------------------------
COPY . /opt/ecolityper/

# --- Editable install of the project itself ----------------------------------
# --no-deps: environment.yml is the single source of truth; pip must not
#            re-resolve bioconda-pinned dependencies.
# -e:       install points at /opt/ecolityper/, so any .py edit takes effect
#            immediately inside the running container (useful for debugging).
RUN micromamba run -n ecolityper pip install --no-deps -e /opt/ecolityper

# --- Use the env as the default shell for subsequent RUN commands ------------
SHELL ["micromamba", "run", "-n", "ecolityper", "/bin/bash", "-c"]

# --- One-time database setup --------------------------------------------------
# abricate: all bundled resistance/virulence/plasmid databases
RUN abricate --setupdb

# AMRFinderPlus: latest reference DB.
# The orchestrator's `--update-amr-db` calls ecoli_amrfinder.py --update-db
# with cwd = the module directory, so we cd there too and invoke directly.
RUN cd /opt/ecolityper/ecoliTyper/modules/amrfinder_module && \
    python ecoli_amrfinder.py --update-db

# E. coli MLST scheme (Achtman + legacy). Baked in so `-u $(id -u):$(id -g)`
# can use MLST without needing to write to the DB directory at runtime.
# This line needs network access at build time — remove it and pull manually
# if building offline.
RUN ecolityper --pull-mlst-db

# --- Runtime environment ------------------------------------------------------
# HOME=/tmp is world-writable so `-u $(id -u):$(id -g)` can write to $HOME.
# Without this, matplotlib, plotly, and amrfinder all crash on first cache write.
ENV HOME=/tmp
ENV MAMBA_ROOT_PREFIX=/opt/conda

# --- Inline entrypoint --------------------------------------------------------
# Sets up the full conda env activation (PATH, LD_LIBRARY_PATH, CONDA_PREFIX)
# that bioconda binaries like amrfinder, diamond, blastn, and prodigal need
# at runtime. Also prepends `ecolityper` if the first argument starts with a
# dash, so `docker run image -i x.fna -o out` works directly.
RUN printf '%s\n' \
    '#!/bin/bash' \
    'set -e' \
    'export CONDA_PREFIX=/opt/conda/envs/ecolityper' \
    'export CONDA_DEFAULT_ENV=ecolityper' \
    'export PATH=${CONDA_PREFIX}/bin:${PATH}' \
    'export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}' \
    'if [ $# -gt 0 ] && [ "${1#-}" != "$1" ]; then' \
    '    set -- ecolityper "$@"' \
    'fi' \
    'exec "$@"' \
    > /usr/local/bin/entrypoint.sh \
    && chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["ecolityper", "-h"]
