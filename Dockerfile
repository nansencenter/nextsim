ARG BASE_IMAGE=nansencenter/nextsim_base:latest
FROM $BASE_IMAGE

# Prepare environment
ENV NEXTSIM_MESH_DIR=/mesh \
    NEXTSIM_DATA_DIR=/data \
    NEXTSIMDIR=/nextsim \
# Compiler names
    CC=mpicc \
    CXX=mpicxx \
    FC=mpifort \
# Compiler options
    MPI_INC_DIR=/usr/lib/openmpi/include \
    CFLAGS="-O3 -fPIC" \
    CXXFLAGS="-O3 -pthread -fPIC -fopenmp "

# Compiler options
ENV CCFLAGS=$CFLAGS

# copy source, compile and copy libs of mapx and bamg
COPY contrib $NEXTSIMDIR/contrib
WORKDIR $NEXTSIMDIR/contrib/mapx/src
RUN make -j8 \
&&  mkdir -p /opt/local/mapx/lib  \
&&  cp -d $NEXTSIMDIR/lib/libmapx* /opt/local/mapx/lib/ \
&&  cp -r $NEXTSIMDIR/contrib/mapx/include /opt/local/mapx \
&&  echo /opt/local/mapx/lib/ >> /etc/ld.so.conf
WORKDIR $NEXTSIMDIR/contrib/bamg/src
RUN make -j8 \
&&  mkdir -p /opt/local/bamg/lib  \
&&  cp -d $NEXTSIMDIR/lib/libbamg* /opt/local/bamg/lib/ \
&&  cp -r $NEXTSIMDIR/contrib/bamg/include /opt/local/bamg \
&&  echo /opt/local/bamg/lib/ >> /etc/ld.so.conf

# copy source, compile and copy libs of core
COPY core $NEXTSIMDIR/core
WORKDIR $NEXTSIMDIR/core/src
RUN make -j8 \
&&  mkdir -p /opt/local/nextsim/lib \
&&  cp -d $NEXTSIMDIR/lib/libnextsim* /opt/local/nextsim/lib \
&&  echo /opt/local/nextsim/lib >> /etc/ld.so.conf

RUN ldconfig

# copy source, compile and copy exec of the model
COPY model $NEXTSIMDIR/model
COPY .git $NEXTSIMDIR/.git
WORKDIR $NEXTSIMDIR/model
RUN make -j8 \
&& cp $NEXTSIMDIR/model/bin/nextsim.exec /usr/local/bin/

RUN rm -rf $NEXTSIMDIR

WORKDIR /root

# allow model to compile in-place
ENV LIBRARY_PATH=/opt/local/bamg/lib:/opt/local/mapx/lib:/opt/local/nextsim/lib \
    PATH=$NEXTSIMDIR/model/bin:$PATH
