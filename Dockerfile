FROM python:3.9.2-buster

MAINTAINER Ai Okada <aokada@ncc.go.jp>

RUN apt-get -y update && \
    apt-get install -y git

WORKDIR /tools

RUN git clone https://github.com/ncc-ccat-gap/GCATWorkflow.git && \
    cd GCATWorkflow && \
    python setup.py build install
