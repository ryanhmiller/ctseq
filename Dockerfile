FROM continuumio/miniconda3
RUN conda create -n ctseqEnv python=3.7
RUN echo "source activate ctseqEnv" > ~/.bashrc

RUN conda install -c bioconda cutadapt
RUN conda install -c bioconda bismark
RUN conda install -c bioconda umi_tools
RUN conda install -c bioconda simplesam

RUN apt-get install git

RUN git clone https://github.com/ryanhmiller/ctseq.git

WORKDIR /ctseq
RUN python setup.py install

WORKDIR /

RUN rm -r ctseq
