FROM ubuntu:latest

RUN apt update
RUN apt-get -y install build-essential
RUN apt-get -y install wget tar gzip

RUN wget http://eddylab.org/software/hmmer/hmmer.tar.gz
RUN tar zxf hmmer.tar.gz
WORKDIR /hmmer-3.2.1
RUN ./configure
RUN make
RUN make install

WORKDIR /conda_dir
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
RUN chmod 0755 miniconda.sh
RUN ./miniconda.sh -b -p /conda_dir/miniconda
ENV PATH="/conda_dir/miniconda/bin:$PATH"
RUN conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda
RUN conda install Biopython
