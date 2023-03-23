FROM python:3.9-bullseye

# Optimise python and debian in docker
ENV PYTHONUNBUFFERED 1
ENV DEBIAN_FRONTEND noninteractive

# Configure /notebooks as a volume
VOLUME /notebooks

# Create the required directories
RUN mkdir -p /app
RUN mkdir -p /vol/web/media
RUN mkdir -p /vol/web/static
RUN mkdir -p /vol/backups

# Update the package cache
RUN apt-get update



# Used for jupyter
COPY ./requirements.jupyter.txt /
ARG ENABLE_JUPYTER
RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then pip install -r /requirements.jupyter.txt; fi
RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then python -m sos_notebook.install; fi
RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then python -m bash_kernel.install; fi
RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then \
cd /tmp && \
wget -q https://github.com/jgm/pandoc/releases/download/2.19.2/pandoc-2.19.2-1-amd64.deb && \
dpkg -i pandoc-2.19.2-1-amd64.deb && \
apt-get update && apt-get install -y \
texlive-xetex \
gconf-service \
libasound2 \
libatk1.0-0 \
libatk-bridge2.0-0 \
libc6 \
libcairo2 \
libcups2 \
libdbus-1-3 \
libexpat1 \
libfontconfig1 \
libgcc1 \
libgconf-2-4 \
libgdk-pixbuf2.0-0 \
libglib2.0-0 \
libgtk-3-0 \
libnspr4 \
libpango-1.0-0 \
libpangocairo-1.0-0 \
libstdc++6 \
libx11-6 \
libx11-xcb1 \
libxcb1 \
libxcomposite1 \
libxcursor1 \
libxdamage1 \
libxext6 \
libxfixes3 \
libxi6 \
libxrandr2 \
libxrender1 \
libxss1 \
libxtst6 \
ca-certificates \
fonts-liberation \
libappindicator1 \
libnss3 \
lsb-release \
xdg-utils \
wget \
nodejs \
npm \
;fi
# RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then jupyter labextension install @techrah/text-shortcuts; fi
RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then jupyter contrib nbextension install; fi

# Jupyter R Kernel
ARG ENABLE_RKERNEL
RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then apt-get install -y r-base; fi
RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then R -e 'install.packages("IRkernel")'; fi
RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then R -e 'IRkernel::installspec()'; fi
RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then R -e 'install.packages("arrow")'; fi
RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then R -e 'install.packages(c("ggplot", "DT", "bioassays"))'; fi

# Used for LDAP authentication
RUN apt-get install -y libsasl2-dev python-dev libldap2-dev libssl-dev

# Global python requirements
COPY ./requirements.txt /
RUN pip --disable-pip-version-check install --upgrade pip
RUN pip install -r /requirements.txt


# Install Node.js and Yarn
RUN curl -sL https://deb.nodesource.com/setup_14.x | bash -
RUN apt-get install -y nodejs
RUN apt-get install jq -y
RUN npm install -g yarn


# Set working directory
WORKDIR /app