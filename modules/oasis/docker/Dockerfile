ARG BASE_IMAGE=nextsim-simple:latest
FROM $BASE_IMAGE

RUN apt-get update --fix-missing && apt-get install -y \
    gfortran \
    libnetcdff-dev \
&&  apt-get clean \
&&  rm -rf /var/lib/apt/lists/* \
&&  ln -s /usr/bin/make /usr/bin/gmake

ENV OASIS_DIR=/opt/local/oasis3-mct \
    USE_OASIS=1
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OASIS_DIR/lib

# get our compilation settings
COPY make.inc make.docker_ubuntu /tmp/

# download OASIS
WORKDIR /tmp
RUN git clone https://gitlab.com/cerfacs/oasis3-mct.git
# compile OASIS
WORKDIR /tmp/oasis3-mct/util/make_dir
RUN cp /tmp/make.* . \
&&  ln -s TopMakefileOasis3 Makefile \
&&  make

WORKDIR $NEXTSIMDIR
RUN rm -r /tmp/* $OASIS_DIR/build
