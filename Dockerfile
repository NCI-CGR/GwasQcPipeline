FROM continuumio/miniconda3
RUN apt-get update

#NOTE: I don't know whether the --fix-broken is actually needed
RUN apt-get -y install gcc libdeflate-dev libcurl4-openssl-dev --fix-broken

WORKDIR /cgr-dev

RUN mkdir GwasQcPipeline

WORKDIR /cgr-dev/GwasQcPipeline

ENV PYTHONPATH=${PYTHONPATH}:${PWD}

COPY . .

#NOTE: I don't think we need a special env inside the container
#RUN conda create -n cgr-dev python=3.8 poetry make -y
#RUN echo "conda activate cgr-dev" >> ~/.bashrc
#RUN conda init bash && conda activate cgr-dev

RUN conda install poetry make -y
RUN conda install -c conda-forge mamba
RUN poetry config virtualenvs.create false
RUN poetry install --without dev

#NOTE: Aren't these also not needed? Unless somehow the docs adds the help for the CLI...
#RUN make -C docs html
#RUN pytest -v
#RUN pre-commit install
#RUN pre-commit run

WORKDIR /cgr-dev

ENTRYPOINT [ "cgr" ]
