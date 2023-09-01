FROM python:3.9

RUN apt-get update

#NOTE: I don't know whether the --fix-broken is actually needed
RUN apt-get -y install gcc make zip libdeflate-dev libcurl4-openssl-dev curl --fix-broken

WORKDIR /home

SHELL ["/bin/bash", "-c"]
RUN wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
RUN bash Mambaforge-$(uname)-$(uname -m).sh -p /home/conda -b
ENV PATH=/home/conda/bin:$PATH
ENV CONDA_EXE=/home/conda/bin/conda
RUN source /home/conda/etc/profile.d/conda.sh && conda activate

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
