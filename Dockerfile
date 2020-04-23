FROM ubuntu:16.04

RUN apt-get update && apt-get install -y \
    git \
    wget \
    unzip \
    build-essential \
    zlib1g-dev 

RUN DEBIAN_FRONTEND=noninteractive apt install -y iproute2 cmake python python-pip build-essential gfortran wget curl
RUN pip install supervisor awscli

ADD manyglucose-4.1-60 manyglucose-4.1-60
RUN cd manyglucose-4.1-60/parallel && make

ADD run.sh supervised-scripts/run.sh
RUN chmod 755 supervised-scripts/run.sh

CMD supervised-scripts/run.sh
