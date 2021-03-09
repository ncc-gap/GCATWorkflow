FROM python:3.9.2-buster

MAINTAINER Ai Okada <aokada@ncc.go.jp>

WORKDIR /tools
RUN apt-get -y update && \
    pip install --upgrade cython && \
    git clone https://github.com/ncc-ccat-gap/GCATWorkflow.git && \
    cd GCATWorkflow && \
    python setup.py build install
