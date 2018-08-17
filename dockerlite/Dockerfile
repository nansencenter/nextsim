FROM boost_petsc_gmsh:boost

RUN apt-get update && apt-get install -y \
    gdb \
    libnetcdf-dev \
    libnetcdf-c++4-dev \
    ssh \
    vim \
    xterm \
    libpetsc3.7.7-dev \
    libgmsh-dev \
    libmumps-ptscotch-dev

ENV NEXTSIMDIR=/nextsim


SHELL ["/bin/bash", "-c"]

COPY Makefile $NEXTSIMDIR/
COPY contrib $NEXTSIMDIR/contrib
COPY core $NEXTSIMDIR/core
WORKDIR $NEXTSIMDIR

ENV BOOST_INCDIR=/opt/local/boost/include
ENV BOOST_LIBDIR=/opt/local/boost/lib

ENV PETSC_DIR=/usr/lib/x86_64-linux-gnu
ENV PETSC_ARCH=arch-linux2-c-opt

ENV GMSH_DIR=/usr/lib/x86_64-linux-gnu

ENV NETCDF_DIR=/usr/lib/x86_64-linux-gnu
ENV OPENMPI_LIB_DIR=/usr/lib/x86_64-linux-gnu/openmpi/lib
ENV OPENMPI_INCLUDE_DIR=/usr/lib/x86_64-linux-gnu/openmpi/include

#RUN make

#COPY model $NEXTSIMDIR/model
#WORKDIR $NEXTSIMDIR/model
#RUN make

#COPY docker/nextsim.sh $NEXTSIMDIR/model/bin/

#RUN ln -s /opt/local/gmsh/bin/gmsh /usr/local/bin/gmsh

#WORKDIR /nextsimdir/docker
#ENTRYPOINT ["/nextsim/model/bin/nextsim.sh"]
CMD [ "/bin/bash" ]
