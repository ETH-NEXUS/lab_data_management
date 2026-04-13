FROM ghcr.io/astral-sh/uv:0.9.10-python3.11-bookworm

# Optimise python and debian in docker
ENV PYTHONUNBUFFERED 1
ENV PIP_BREAK_SYSTEM_PACKAGES 1
ENV PIP_DISABLE_PIP_VERSION_CHECK 1
ENV DEBIAN_FRONTEND noninteractive
ENV UV_COMPILE_BYTECODE 1
ENV UV_LINK_MODE copy

# Configure /notebooks as a volume
VOLUME /notebooks

# Create the required directories
RUN mkdir -p /app
RUN mkdir -p /vol/web/media
RUN mkdir -p /vol/web/static
RUN mkdir -p /vol/backups

# Update the package cache
RUN apt-get update --fix-missing

# Install PostgreSQL client
RUN apt-get install -y postgresql-client git

# Cleanup
RUN apt-get clean && rm -rf /var/lib/apt/lists/*


# Install Node.js, npm and Yarn
# According to this guide: 
# https://github.com/nodesource/distributions
# RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash -
# RUN apt-get install -y nodejs jq
# RUN npm install -g yarn
# Since this is deprecated we use the newer way according to the following information
# https://github.com/nodesource/distributions#debian-and-ubuntu-based-distributions
RUN curl -fsSL https://deb.nodesource.com/gpgkey/nodesource-repo.gpg.key | gpg --dearmor -o /etc/apt/keyrings/nodesource.gpg
RUN echo "deb [signed-by=/etc/apt/keyrings/nodesource.gpg] https://deb.nodesource.com/node_20.x nodistro main" | tee /etc/apt/sources.list.d/nodesource.list
RUN apt-get update --fix-missing && apt-get install nodejs -y

# Install R
# According to this guide: 
# https://www.linuxcapable.com/how-to-install-r-programming-language-on-debian-linux/
ARG ENABLE_R
RUN if [ "${ENABLE_R}" = "True" ]; then \
apt-get install -y dirmngr apt-transport-https ca-certificates software-properties-common && \
gpg --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' && \
gpg --armor --export '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' | gpg --dearmor | tee /usr/share/keyrings/cran.gpg > /dev/null && \
echo "deb [signed-by=/usr/share/keyrings/cran.gpg] https://cloud.r-project.org/bin/linux/debian bookworm-cran40/" | tee /etc/apt/sources.list.d/cran.list \
;fi
RUN if [ "${ENABLE_R}" = "True" ]; then apt-get update && apt-get install -y r-base cmake; fi
# R dependencies are installed at the end to avoid rebuild the whole image in case we want to add 
# an additional R dependency

ARG ENABLE_JUPYTER
ARG TARGETARCH

RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then \
      ARCH="${TARGETARCH:-$(dpkg --print-architecture)}"; \
      case "${ARCH}" in \
        amd64) \
          QUARTO_DEB="quarto-1.4.156-linux-amd64.deb"; \
          wget -q "https://github.com/quarto-dev/quarto-cli/releases/download/v1.4.156/${QUARTO_DEB}"; \
          dpkg -i "${QUARTO_DEB}"; \
          quarto install tinytex; \
          ;; \
        arm64) \
          QUARTO_DEB="quarto-1.4.156-linux-arm64.deb"; \
          wget -q "https://github.com/quarto-dev/quarto-cli/releases/download/v1.4.156/${QUARTO_DEB}"; \
          dpkg -i "${QUARTO_DEB}"; \
          wget -qO- "https://tinytex.yihui.org/install-bin-unix.sh" | sh; \
          ;; \
        *) \
          echo "Unsupported architecture for Quarto: ${ARCH}" >&2; \
          exit 1; \
          ;; \
      esac; \
    fi
ENV PATH="/root/bin:/root/.TinyTeX/bin/aarch64-linux:/root/.TinyTeX/bin/x86_64-linux:${PATH}"

# Used for LDAP authentication
RUN apt-get install -y libsasl2-dev python-dev-is-python3 libldap2-dev libssl-dev



# 2. Optimized R Layer
RUN if [ "${ENABLE_R}" = "True" ]; then \
      R -q -e 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/latest")); install.packages(c( \
        "ggplot2", "DT", "BioRssay", "repr", "tidyverse", "car", \
        "ggpubr", "lawstat", "RColorBrewer", "gplots", "ggrepel", \
        "cowplot", "ggExtra", "pheatmap", "reshape2", "drc", \
        "knitr", "locfit", "robust", "lubridate", "rmarkdown" \
      ), Ncpus = max(1L, parallel::detectCores() - 1L))'; \
    fi

# RUN --mount=type=cache,target=/root/.cache/pak \
#     if [ "${ENABLE_R}" = "True" ]; then \
#       R -q -e 'install.packages("pak", repos = "https://r-lib.github.io/p/pak/stable")'; \
#       R -q -e ' \
#         options(pak.no_extra_messages = TRUE); \
#         options(pak.num_workers = parallel::detectCores()); \
#         pkgs <- c("ggplot2", "DT", "BioRssay", "repr", "tidyverse", "car", "ggpubr", "lawstat", "RColorBrewer", "gplots", "ggrepel", "cowplot", "ggExtra", "pheatmap", "reshape2", "drc", "knitr", "locfit", "robust", "lubridate", "rmarkdown"); \
#         sysreqs <- pak::pkg_sysreqs(pkgs); \
#         system(paste(c("apt-get update -qq", sysreqs$install_scripts), collapse = " && ")); \
#         pak::pkg_install(pkgs); \
#       '; \
#     fi


COPY ./pyproject.toml ./uv.lock /
ARG DEPENDENCIES=prod
RUN --mount=type=cache,target=/root/.cache/uv \
    /bin/sh -eux -c '\
      UV_ARGS="--frozen --no-install-project --no-default-groups"; \
      if [ "$DEPENDENCIES" = "dev" ]; then \
        UV_ARGS="$UV_ARGS --group dev"; \
      fi; \
      if [ "$ENABLE_JUPYTER" = "True" ]; then \
        UV_ARGS="$UV_ARGS --group jupyter"; \
      fi; \
      if [ "$ENABLE_R" = "True" ]; then \
        UV_ARGS="$UV_ARGS --group r"; \
      fi; \
      uv sync $UV_ARGS; \
    '
ENV PATH="/.venv/bin:$PATH"

# Set working directory
WORKDIR /app
