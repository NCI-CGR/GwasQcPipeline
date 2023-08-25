#FROM continuumio/miniconda3:4.11.0
FROM condaforge/mambaforge:4.11.0-0

RUN apt-get update

#NOTE: I don't know whether the --fix-broken is actually needed
RUN apt-get -y install gcc make zip libdeflate-dev libcurl4-openssl-dev curl --fix-broken

WORKDIR /home

RUN mkdir GwasQcPipeline
RUN mkdir data

WORKDIR /home/GwasQcPipeline

ENV PYTHONPATH=${PYTHONPATH}:${PWD}

COPY . .

RUN curl -sSL https://install.python-poetry.org | python3
ENV PATH=/root/.local/bin:${PATH}
RUN poetry config virtualenvs.create false
RUN poetry install --without dev

WORKDIR /home/data

#ENTRYPOINT [ "cgr" ]
ENTRYPOINT [ "/bin/bash" ]
