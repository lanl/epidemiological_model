FROM continuumio/miniconda3

# make sure the shell is always a bash login shell
SHELL ["/bin/bash", "--login", "-c"]

# install build essentials, which will include things like gcc (needed to build the R smooth package below)
RUN apt-get update

# configure the proxy
RUN echo "export HTTP_PROXY=http://proxyout.lanl.gov:8080" >> ~/.bashrc
RUN echo "export http_proxy=$HTTP_PROXY" >> ~/.bashrc
RUN echo "export HTTPS_PROXY=$HTTP_PROXY" >> ~/.bashrc
RUN echo "export https_proxy=$HTTP_PROXY" >> ~/.bashrc
RUN echo "export NO_PROXY=localhost" >> ~/.bashrc

# setup the conda environment
RUN conda init bash
RUN conda create -y -n human-epi-env
RUN echo "conda activate human-epi-env" >> ~/.bashrc
RUN conda install -y -c conda-forge numpy pyyaml pandas scipy pyarrow matplotlib sphinx pytest
