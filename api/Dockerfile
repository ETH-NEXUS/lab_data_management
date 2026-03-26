FROM ghcr.io/astral-sh/uv:0.9.10-python3.13-bookworm AS api-deps

ARG DEPENDENCIES=prod
ARG ENABLE_JUPYTER
ARG ENABLE_R

ENV PYTHONUNBUFFERED=1 \
    PIP_BREAK_SYSTEM_PACKAGES=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    DEBIAN_FRONTEND=noninteractive \
    UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy

VOLUME /notebooks

# Create required directories
RUN mkdir -p /app /vol/web/media /vol/web/static /vol/backups

# System dependencies
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
      postgresql-client \
      libsasl2-dev \
      python-dev-is-python3 \
      libldap2-dev \
      libssl-dev \
      wget \
      gnupg \
      dirmngr \
      apt-transport-https \
      ca-certificates \
      software-properties-common \
      cmake \
      git && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install latest available R for Debian Bookworm from CRAN
RUN if [ "${ENABLE_R}" = "True" ]; then \
      gpg --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' && \
      gpg --armor --export '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' | gpg --dearmor > /usr/share/keyrings/cran.gpg && \
      echo "deb [signed-by=/usr/share/keyrings/cran.gpg] https://cloud.r-project.org/bin/linux/debian bookworm-cran40/" > /etc/apt/sources.list.d/cran.list && \
      apt-get update && \
      apt-get install -y --no-install-recommends r-base && \
      apt-get clean && rm -rf /var/lib/apt/lists/*; \
    fi

# Jupyter tooling
RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then \
      wget -q https://github.com/quarto-dev/quarto-cli/releases/download/v1.4.156/quarto-1.4.156-linux-amd64.deb && \
      dpkg -i quarto-1.4.156-linux-amd64.deb && \
      rm -f quarto-1.4.156-linux-amd64.deb && \
      quarto install tinytex; \
    fi

# Install R dependencies
RUN if [ "${ENABLE_R}" = "True" ]; then R -e 'install.packages(c("ggplot2", "DT", "bioassays", "repr", "tidyverse"), dependencies = TRUE, INSTALL_opts = c("--no-lock"))'; fi
RUN if [ "${ENABLE_R}" = "True" ]; then R -e 'install.packages(c("car", "ggpubr", "captioner", "lawstat"), dependencies = TRUE, INSTALL_opts = c("--no-lock"))'; fi
RUN if [ "${ENABLE_R}" = "True" ]; then R -e 'install.packages(c("RColorBrewer", "gplots", "ggrepel", "cowplot", "ggExtra", "pheatmap", "reshape2", "drc", "knitr", "locfit", "robust"), dependencies = TRUE, INSTALL_opts = c("--no-lock"))'; fi
RUN if [ "${ENABLE_R}" = "True" ]; then R -e 'install.packages(c("lubridate", "rmarkdown"), dependencies = TRUE, INSTALL_opts = c("--no-lock"))'; fi

WORKDIR /

# Install python dependencies with uv
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    /bin/sh -eux -c '\
      UV_ARGS="--no-install-project --no-default-groups"; \
      if [ "$DEPENDENCIES" = "dev" ]; then \
        UV_ARGS="$UV_ARGS --group dev --group jupyter"; \
      else \
        if [ "$ENABLE_JUPYTER" = "True" ]; then \
          UV_ARGS="$UV_ARGS --group jupyter"; \
        fi; \
      fi; \
      if [ "$ENABLE_R" = "True" ]; then \
        UV_ARGS="$UV_ARGS --group r"; \
      fi; \
      uv sync $UV_ARGS; \
    '

ENV PATH="/.venv/bin:$PATH"

FROM api-deps AS api-runtime

COPY ./app /app
RUN chmod +x /app/entrypoint.sh

WORKDIR /app
