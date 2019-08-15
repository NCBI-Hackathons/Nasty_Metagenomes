FROM ubuntu:latest

RUN apt update
RUN apt-get -y install build-essential
RUN apt-get -y install wget tar gzip

RUN wget http://eddylab.org/software/hmmer/hmmer.tar.gz
RUN tar zxf hmmer.tar.gz
WORKDIR hmmer-3.2.1
RUN ./configure
RUN make
RUN make install 
