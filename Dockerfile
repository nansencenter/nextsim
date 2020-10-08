ARG BASE_IMAGE=nansencenter/nextsim_base_dev:latest
FROM $BASE_IMAGE

ENV OPENMPI_INCLUDE_DIR=/usr/lib/x86_64-linux-gnu/openmpi/include \
    OPENMPI_LIB_DIR=/usr/lib/x86_64-linux-gnu/openmpi/lib \
    BOOST_INCDIR=/opt/local/boost/include \
    BOOST_LIBDIR=/opt/local/boost/lib \
    PETSC_DIR=/opt/local/petsc \
    GMSH_DIR=/opt/local/gmsh \
    NEXTSIM_MESH_DIR=/mesh \
    NEXTSIM_DATA_DIR=/data \
    LIBRARY_PATH=/opt/local/mapx/lib:/opt/local/bamg/lib:/opt/local/nextsim/lib

# copy source code into temporary dir
ENV NEXTSIMDIR=/tmp/nextsim
COPY contrib $NEXTSIMDIR/contrib
COPY core $NEXTSIMDIR/core
COPY model $NEXTSIMDIR/model
COPY modules $NEXTSIMDIR/modules
COPY Makefile $NEXTSIMDIR/Makefile

# compile mapx
WORKDIR $NEXTSIMDIR/contrib/mapx/src
RUN make -j8 \
&&  mkdir -p /opt/local/mapx/lib  \
&&  cp -d $NEXTSIMDIR/lib/libmapx* /opt/local/mapx/lib/ \
&&  cp -r $NEXTSIMDIR/contrib/mapx/include /opt/local/mapx \
&&  echo /opt/local/mapx/lib/ >> /etc/ld.so.conf \
&&  ldconfig

# compile bamg
WORKDIR $NEXTSIMDIR/contrib/bamg/src
RUN make -j8 \
&&  mkdir -p /opt/local/bamg/lib  \
&&  cp -d $NEXTSIMDIR/lib/libbamg* /opt/local/bamg/lib/ \
&&  cp -r $NEXTSIMDIR/contrib/bamg/include /opt/local/bamg \
&&  echo /opt/local/bamg/lib/ >> /etc/ld.so.conf \
&&  ldconfig

# compile model
WORKDIR $NEXTSIMDIR
RUN make docker && \
    mkdir -p /opt/local/nextsim/lib && \
    cp -d $NEXTSIMDIR/lib/* /opt/local/nextsim/lib/ && \
    echo /opt/local/nextsim/lib/ >> /etc/ld.so.conf && \
    ldconfig && \
    cp $NEXTSIMDIR/model/bin/nextsim.exec /usr/local/bin/ && \
    rm -rf $NEXTSIMDIR

# trick to use inplace compiled nextsim
ENV NEXTSIMDIR=/nextsim
ENV PATH=$NEXTSIMDIR/model/bin:$PATH
WORKDIR $NEXTSIMDIR
