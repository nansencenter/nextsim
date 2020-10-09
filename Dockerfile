ARG BASE_IMAGE=nansencenter/nextsim_base_dev:latest
FROM $BASE_IMAGE

ENV NEXTSIM_MESH_DIR=/mesh \
    NEXTSIM_DATA_DIR=/data \
    NEXTSIMDIR=/nextsim \
    LD_LIBRARY_PATH=/nextsim/lib/ \
    PATH=/nextsim/model/bin:$PATH

# copy source code into temporary dir
COPY contrib $NEXTSIMDIR/contrib
COPY core $NEXTSIMDIR/core
COPY model $NEXTSIMDIR/model
COPY modules $NEXTSIMDIR/modules
COPY Makefile $NEXTSIMDIR/Makefile

# compile model
WORKDIR $NEXTSIMDIR
RUN make -j8

# copy mapx
WORKDIR $NEXTSIMDIR/contrib/mapx/src
RUN mkdir -p /opt/local/mapx/lib  \
&&  cp -d $NEXTSIMDIR/lib/libmapx* /opt/local/mapx/lib/ \
&&  cp -r $NEXTSIMDIR/contrib/mapx/include /opt/local/mapx

# copy bamg
WORKDIR $NEXTSIMDIR/contrib/bamg/src
RUN mkdir -p /opt/local/bamg/lib  \
&&  cp -d $NEXTSIMDIR/lib/libbamg* /opt/local/bamg/lib/ \
&&  cp -r $NEXTSIMDIR/contrib/bamg/include /opt/local/bamg

WORKDIR /root
