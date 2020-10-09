ARG BASE_IMAGE=nansencenter/nextsim_base_dev:latest
FROM $BASE_IMAGE

ENV NEXTSIM_MESH_DIR=/mesh \
    NEXTSIM_DATA_DIR=/data \
    LIBRARY_PATH=/opt/local/mapx/lib:/opt/local/bamg/lib \
    NEXTSIMDIR=/nextsim \
    LD_LIBRARY_PATH=/nextsim/lib/ \
    PATH=/nextsim/model/bin:$PATH

# copy source code into temporary dir
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
RUN make docker
