FROM postgres:15-bullseye

RUN apt-get update && apt-get install -y postgresql-plpython3-15 python3-pip

# Global python requirements
COPY ./requirements.txt /
RUN pip --disable-pip-version-check install --upgrade pip
RUN pip install -r /requirements.txt