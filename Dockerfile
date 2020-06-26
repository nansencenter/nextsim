FROM nansencenter/boost_petsc_gmsh:0.0.6

RUN apt-get update && apt-get install -yq --no-install-recommends \
    libnetcdf-dev \
    libnetcdf-c++4-dev \
    ssh \
    valgrind \
    bc \
    git \
&&  apt-get clean \
&&  rm -rf /var/lib/apt/lists/*

ENV OPENMPI_INCLUDE_DIR=/usr/lib/x86_64-linux-gnu/openmpi/include \
    OPENMPI_LIB_DIR=/usr/lib/x86_64-linux-gnu/openmpi/lib \
    BOOST_INCDIR=/opt/local/boost/include \
    BOOST_LIBDIR=/opt/local/boost/lib \
    PETSC_DIR=/opt/local/petsc \
    GMSH_DIR=/opt/local/gmsh \
    PATH=$PATH:/nextsim/model/bin \
    NEXTSIMDIR=/nextsim \
    NEXTSIM_MESH_DIR=/mesh \
    NEXTSIM_DATA_DIR=/data \
    LIBRARY_PATH=/opt/local/mapx/lib:/opt/local/bamg/lib

# copy mapx source, compile and copy libs and include into
# /opt/local/mapx/lib and /opt/local/mapx/include
COPY contrib/mapx $NEXTSIMDIR/contrib/mapx
WORKDIR $NEXTSIMDIR/contrib/mapx/src

RUN make -j8 \
&&  mkdir -p /opt/local/mapx/lib  \
&&  cp -d $NEXTSIMDIR/lib/libmapx* /opt/local/mapx/lib/ \
&&  cp -r $NEXTSIMDIR/contrib/mapx/include /opt/local/mapx \
&&  echo /opt/local/mapx/lib/ >> /etc/ld.so.conf \
&&  ldconfig \
&&  rm -rf $NEXTSIMDIR/contrib/mapx

# copy bamg source, compile and copy libs and include into
# /opt/local/bamg/lib and /opt/local/bamg/include
COPY contrib/bamg $NEXTSIMDIR/contrib/bamg
WORKDIR $NEXTSIMDIR/contrib/bamg/src
RUN make -j8 \
&&  mkdir -p /opt/local/bamg/lib  \
&&  cp -d $NEXTSIMDIR/lib/libbamg* /opt/local/bamg/lib/ \
&&  cp -r $NEXTSIMDIR/contrib/bamg/include /opt/local/bamg \
&&  echo /opt/local/bamg/lib/ >> /etc/ld.so.conf \
&&  ldconfig \
&&  rm -rf $NEXTSIMDIR/contrib/bamg

WORKDIR $NEXTSIMDIR
