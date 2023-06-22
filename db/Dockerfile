FROM postgres:15-bookworm

ENV PIP_BREAK_SYSTEM_PACKAGES 1
ENV PIP_DISABLE_PIP_VERSION_CHECK 1

RUN apt-get update && apt-get install -y postgresql-plpython3-15 python3-pip

# Global python requirements
COPY ./requirements.txt /
RUN pip install --upgrade pip
RUN pip install -r /requirements.txt