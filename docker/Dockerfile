FROM continuumio/miniconda3:main

ADD env.yml .

RUN conda env create -f env.yml
RUN echo "source activate RNASeq" > ~/.bashrc

ENV PATH /opt/conda/envs/RNASeq/bin:$PATH
