FROM akorosov/boost_petsc_gmsh:0.0.2

RUN apt-get update && apt-get install -y \
    libnetcdf-dev \
    libnetcdf-c++4-dev \
    ssh \
&& rm -rf /var/lib/apt/lists/*

ENV NEXTSIMDIR=/nextsim


SHELL ["/bin/bash", "-c"]

COPY Makefile $NEXTSIMDIR/
COPY contrib $NEXTSIMDIR/contrib
COPY core $NEXTSIMDIR/core
WORKDIR $NEXTSIMDIR
RUN source $NEXTSIMDIR/nextsim.src && make

COPY model $NEXTSIMDIR/model
WORKDIR $NEXTSIMDIR/model
RUN source $NEXTSIMDIR/nextsim.src && make

COPY docker/nextsim.sh $NEXTSIMDIR/model/bin/
COPY docker/test.cfg $NEXTSIMDIR/model/bin/

ENTRYPOINT ["/nextsim/model/bin/nextsim.sh"]
#CMD [ "/bin/bash" ]
