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

# Jupyter R Kernel
ARG ENABLE_RKERNEL
RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then apt-get install -y dirmngr gnupg apt-transport-https ca-certificates software-properties-common; fi
# The foolowing two lines is to upgrade R to the latest version on bullseye
# Reference: https://linux.how2shout.com/install-the-latest-r-programming-language-version-on-debian-11-bullseye/
RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then apt-key adv --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'; fi
RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then add-apt-repository 'deb http://cloud.r-project.org/bin/linux/debian bullseye-cran40/'; fi
RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then apt-get update && apt-get install -y cmake r-base; fi
# The following packages are needed the Jupyter Kernel to run
RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then R -e 'install.packages("IRkernel")'; fi
RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then R -e 'IRkernel::installspec()'; fi
RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then R -e 'install.packages("arrow")'; fi
# Below the packages that are needed can be configured
RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then R -e 'install.packages(c("ggplot", "DT", "bioassays"), dependencies = TRUE, INSTALL_opts = c("--no-lock"))'; fi
RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then R -e 'install.packages(c("car", "ggpubr"), dependencies = TRUE, INSTALL_opts = c("--no-lock"))'; fi

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
wget -q https://github.com/quarto-dev/quarto-cli/releases/download/v1.3.340/quarto-1.3.340-linux-amd64.deb && \
dpkg -i quarto-1.3.340-linux-amd64.deb && \
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
RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then jupyter contrib nbextension install --user; fi
RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then jupyter nbextension enable printview/main ; fi
RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then jupyter nbextension enable hide_input_all/main ; fi
RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then jupyter nbextension enable hide_input/main ; fi
RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then jupyter nbextension enable livemdpreview/livemdpreview ; fi
RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then jupyter nbextension enable code_prettify/autopep8 ; fi
RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then jupyter nbextension enable varInspector/main ; fi

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