FROM python:3.11-bookworm

# Optimise python and debian in docker
ENV PYTHONUNBUFFERED 1
ENV PIP_BREAK_SYSTEM_PACKAGES 1
ENV PIP_DISABLE_PIP_VERSION_CHECK 1
ENV DEBIAN_FRONTEND noninteractive

# Configure /notebooks as a volume
VOLUME /notebooks

# Create the required directories
RUN mkdir -p /app
RUN mkdir -p /vol/web/media
RUN mkdir -p /vol/web/static
RUN mkdir -p /vol/backups

# Update the package cache
RUN apt-get update --fix-missing

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
# According to this guide: 
# https://anderfernandez.com/en/blog/how-to-program-with-python-and-r-in-the-same-jupyter-notebook/
RUN if [ "${ENABLE_R}" = "True" ]; then pip install rpy2[all]==3.5.12; fi
# R dependencies are installed at the end to avoid rebuild the whole image in case we want to add 
# an additional R dependency

# Update pip
RUN pip install --upgrade pip

# Used for Jupyter
COPY ./requirements.jupyter.txt /
ARG ENABLE_JUPYTER
RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then pip install -r /requirements.jupyter.txt; fi
# RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then python -m sos_notebook.install; fi
# RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then python -m bash_kernel.install; fi

# Build jupyter lab to integrate installed extensions
# RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then jupyter lab build; fi

# Debian dependencies for Jupyter
# RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then \
# cd /tmp && \
# apt-get update && apt-get install -y \
# texlive-xetex \
# gconf-service \
# libasound2 \
# libatk1.0-0 \
# libatk-bridge2.0-0 \
# libc6 \
# libcairo2 \
# libcups2 \
# libdbus-1-3 \
# libexpat1 \
# libfontconfig1 \
# libgcc1 \
# libgconf-2-4 \
# libgdk-pixbuf2.0-0 \
# libglib2.0-0 \
# libgtk-3-0 \
# libnspr4 \
# libpango-1.0-0 \
# libpangocairo-1.0-0 \
# libstdc++6 \
# libx11-6 \
# libx11-xcb1 \
# libxcb1 \
# libxcomposite1 \
# libxcursor1 \
# libxdamage1 \
# libxext6 \
# libxfixes3 \
# libxi6 \
# libxrandr2 \
# libxrender1 \
# libxss1 \
# libxtst6 \
# ca-certificates \
# fonts-liberation \
# libappindicator1 \
# libnss3 \
# lsb-release \
# xdg-utils \
# wget \
# ;fi
# RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then \
# wget -q https://github.com/jgm/pandoc/releases/download/3.1.3/pandoc-3.1.3-1-amd64.deb && \
# dpkg -i pandoc-3.1.3-1-amd64.deb && \
# ;fi
RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then \
wget -q https://github.com/quarto-dev/quarto-cli/releases/download/v1.4.156/quarto-1.4.156-linux-amd64.deb && \
dpkg -i quarto-1.4.156-linux-amd64.deb && \
quarto install tinytex \
;fi
# RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then \
# apt-get install -y xfonts-75dpi xfonts-base && \
# wget -q https://github.com/wkhtmltopdf/packaging/releases/download/0.12.6.1-3/wkhtmltox_0.12.6.1-3.bookworm_amd64.deb && \
# dpkg -i wkhtmltox_0.12.6.1-3.bookworm_amd64.deb \
# ;fi
# RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then pip install hide_code[lab]; fi
# RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then jupyter labextension install @techrah/text-shortcuts; fi
# RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then jupyter contrib nbextension install --user; fi
# RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then jupyter nbextension enable printview/main ; fi
# RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then jupyter nbextension enable hide_input_all/main ; fi
# RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then jupyter nbextension enable hide_input/main ; fi
# RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then jupyter nbextension enable livemdpreview/livemdpreview ; fi
# # RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then jupyter nbextension enable code_prettify/autopep8 ; fi
# RUN if [ "${ENABLE_JUPYTER}" = "True" ]; then jupyter nbextension enable varInspector/main ; fi


# Jupyter R Kernel
# ARG ENABLE_RKERNEL
# RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then apt-get install -y dirmngr gnupg apt-transport-https ca-certificates software-properties-common; fi
# # The foolowing two lines is to upgrade R to the latest version on bullseye
# # Reference: https://linux.how2shout.com/install-the-latest-r-programming-language-version-on-debian-11-bullseye/
# RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then apt-key adv --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'; fi
# RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then add-apt-repository 'deb http://cloud.r-project.org/bin/linux/debian bullseye-cran40/'; fi
# RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then apt-get update && apt-get install -y cmake r-base; fi
# # The following packages are needed the Jupyter Kernel to run
# RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then R -e 'install.packages("IRkernel")'; fi
# RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then R -e 'IRkernel::installspec()'; fi
# RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then R -e 'install.packages("arrow")'; fi
# # Below the packages that are needed can be configured
# RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then R -e 'install.packages(c("ggplot", "DT", "bioassays"), dependencies = TRUE, INSTALL_opts = c("--no-lock"))'; fi
# RUN if [ "${ENABLE_RKERNEL}" = "True" ]; then R -e 'install.packages(c("car", "ggpubr"), dependencies = TRUE, INSTALL_opts = c("--no-lock"))'; fi

# Used for LDAP authentication
RUN apt-get install -y libsasl2-dev python-dev-is-python3 libldap2-dev libssl-dev

# Global python requirements
COPY ./requirements.txt /
RUN pip install -r /requirements.txt

# Install R dependecies
RUN if [ "${ENABLE_R}" = "True" ]; then R -e 'install.packages(c("ggplot2", "DT", "bioassays"), dependencies = TRUE, INSTALL_opts = c("--no-lock"))'; fi
RUN if [ "${ENABLE_R}" = "True" ]; then R -e 'install.packages(c("car", "ggpubr", "captioner", "lawstat"), dependencies = TRUE, INSTALL_opts = c("--no-lock"))'; fi

# Set working directory
WORKDIR /app