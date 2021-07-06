FROM ubuntu:latest
MAINTAINER "Tamas K Stenczel <tks@polychord.co.uk>"

# Time Zone data
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/London

RUN apt-get -q update && \
    apt-get -qy install --no-install-recommends \
        # general
        make \
        gfortran \
        libopenmpi-dev \
        libblas-dev \
        liblapack-dev \
        libfftw3-3  \
        # cmake and tools
        cmake vim git wget ca-certificates \
        g++ \
    && rm -rf /var/lib/apt/lists/*


# DEV use: oh-my-zsh instead of bash
# Uses "robbyrussell" theme (original Oh My Zsh theme), with no plugins
RUN sh -c "$(wget -O- https://github.com/deluan/zsh-in-docker/releases/download/v1.1.1/zsh-in-docker.sh)" -- \
    -t robbyrussell
ENTRYPOINT "/usr/bin/zsh"
WORKDIR /work
