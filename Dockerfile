FROM reslp/mamba:0.15.3
MAINTAINER <philipp.resl@uni-graz.at>

RUN apt-get update --allow-releaseinfo-change && apt-get install -y \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    pkg-config

ENV VERSION=1.17.2 
ENV OS=linux 
ENV ARCH=amd64
RUN  wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
    tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \
    rm go$VERSION.$OS-$ARCH.tar.gz
ENV GOPATH=${HOME}/go
ENV PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin
#RUN echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
#    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
#    source ~/.bashrc

ENV VERSION=v3.6.3
RUN  mkdir -p $GOPATH/src/github.com/sylabs && \
    cd $GOPATH/src/github.com/sylabs && \
    git clone https://github.com/sylabs/singularity.git && \
    cd ./singularity && \
    git checkout $VERSION && \
    ./mconfig && make -C ./builddir && make -C ./builddir install

RUN mamba install -y -c bioconda snakemake=6.0.2

RUN git clone --recursive https://github.com/reslp/phylociraptor.git
WORKDIR /phylociraptor
ENV PATH=${PATH}:/phylociraptor

CMD ["phylociraptor"]

