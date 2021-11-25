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
RUN wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
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

# need to be passed as build arguments to set the correct user and group id
RUN find /opt/conda -type d -exec chmod a+w {} \;

ARG USER_ID
ARG GROUP_ID
ARG USER
ARG GROUP
ARG HOSTTYPE

RUN if [ "$HOSTTYPE" = "darwin" ] ; then echo "MacOS host, will skip"; else addgroup --gid $GROUP_ID $USER; fi
#RUN if [ "$HOSTTYPE" = "darwin" ] ; then echo "MacOS host, will skip"; else adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID $USER; fi
RUN if [ "$HOSTTYPE" = "darwin" ] ; then echo "MacOS host, will skip"; else useradd -ms /bin/bash --uid $USER_ID --gid $GROUP_ID $USER; fi
USER $USER

#RUN chown -R $USER:$USER /opt/conda

RUN mamba create -y -c bioconda -c conda-forge -n snakemake snakemake=6.0.2

WORKDIR /home/$USER
RUN git clone --recursive https://github.com/reslp/phylociraptor.git
WORKDIR /home/$USER/phylociraptor
ENV PATH=${PATH}:/home/$USER/phylociraptor

# need to install bc which is not there by default
RUN apt-get install -y bc

SHELL ["conda", "run", "-n", "snakemake", "/bin/bash", "-c"]

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "snakemake"]
CMD ["phylociraptor"]
