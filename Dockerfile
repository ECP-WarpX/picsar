#copied with slight modifications from Warp repository
FROM ubuntu:latest

# Install a few packages, as root
RUN apt-get update \
    && apt-get install -y \
    sudo \
    wget \
    make \
    git \
    gcc \
    nano \
    gfortran \
    libx11-dev \
    openmpi-bin libopenmpi-dev \
    python3 \
    python3-pip \
    python3-numpy \
    python3-scipy \
    python3-mpi4py \
    python3-h5py \
    libfftw3-dev \ 
    libfftw3-mpi-dev \
    libboost-math-dev \ 
    libboost-test-dev \ 
    virtualenv \
    python3-virtualenv \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# Install pygist
ENV GIT_SSL_NO_VERIFY 1
RUN git clone https://bitbucket.org/dpgrote/pygist.git \
    && cd pygist \
    && python3 setup.py config \
    && python3 setup.py install \
    && cd ../ \
    && rm -rf pygist

# Create a new user and clone the current branch of Warp
# into the Docker container
RUN useradd --create-home warp_user
RUN cd /home/warp_user &&  git clone https://bitbucket.org/berkeleylab/warp/src warp
RUN chown -R warp_user /home/warp_user/warp/
# Grant sudo access without password
RUN echo 'warp_user ALL=(ALL) NOPASSWD: ALL' >> /etc/sudoers

# openPMD-viewer is installed mainly for tests
# Note: matplotlib is installed with pip since the apt-get install matplotlib
#       needs the time zone to be set.
SHELL ["/bin/bash", "-c"]
RUN su warp_user && cd /home/warp_user && virtualenv -p python3 mypy && source mypy/bin/activate && \
    pip3 --no-cache-dir install --upgrade pip && \
    pip3 --no-cache-dir install matplotlib \
    openPMD-viewer \
    matplotlib \
    pytest \
    Forthon

# Switch to the new user
WORKDIR /home/warp_user
USER warp_user

# Compile warp
RUN source mypy/bin/activate \
    && cd warp/pywarp90 \
    && echo 'FCOMP= -F gfortran' >> Makefile.local3 \
    && echo 'FCOMP= -F gfortran' >> Makefile.local3.pympi \
    && echo 'FCOMPEXEC= --fcompexec mpifort' >> Makefile.local3.pympi \
    && make install3 INSTALLOPTIONS=--user \
    && make clean3 \
    && make pinstall3 INSTALLOPTIONS=--user \
    && make pclean3

# This is needed to get around a bug in openmpi that would print copious error messages
# Unfortunately, this turns off CMA and uses shared memory for communication.
# An alternative is to do "docker run --cap-add SYS_PTRACE ...", which keeps CMA.
ENV OMPI_MCA_btl_vader_single_copy_mechanism none

# Prepare the run directory
RUN mkdir run/
WORKDIR /home/warp_user/run/

#change ownership of mypy dir
RUN sudo chown -R warp_user /home/warp_user/mypy
