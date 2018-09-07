FROM akorosov/boost_petsc_gmsh:0.0.5

RUN apt-get update && apt-get install -y \
    libnetcdf-dev \
    libnetcdf-c++4-dev \
    ssh \
&& rm -rf /var/lib/apt/lists/*

ENV NEXTSIMDIR=/opt/local/nextsim
SHELL ["/bin/bash", "-c"]

# copy mapx source, compile and copy libs and include into
# /opt/local/mapx/lib and /opt/local/mapx/include
COPY contrib/mapx $NEXTSIMDIR/contrib/mapx
WORKDIR $NEXTSIMDIR/contrib/mapx/src
RUN source /root/.nextsimrc \
&&  make -j8 \
&&  mkdir -p /opt/local/mapx/lib  \
&&  cp $NEXTSIMDIR/lib/libmapx* /opt/local/mapx/lib/ \
&&  cp -r $NEXTSIMDIR/contrib/mapx/include /opt/local/mapx

# copy bamg source, compile and copy libs and include into
# /opt/local/bamg/lib and /opt/local/bamg/include
COPY contrib/bamg $NEXTSIMDIR/contrib/bamg
WORKDIR $NEXTSIMDIR/contrib/bamg/src
RUN source /root/.nextsimrc \
&&  make -j8 \
&&  mkdir -p /opt/local/bamg/lib  \
&&  cp $NEXTSIMDIR/lib/libbamg* /opt/local/bamg/lib/ \
&&  cp -r $NEXTSIMDIR/contrib/bamg/include /opt/local/bamg

# copy core source and compile
COPY core $NEXTSIMDIR/core
WORKDIR $NEXTSIMDIR/core/src
RUN source /root/.nextsimrc && make -j8

# copy model source and compile
COPY model $NEXTSIMDIR/model
WORKDIR $NEXTSIMDIR/model
RUN source /root/.nextsimrc && make -j8

# install nextsim system-wide
RUN ln -s $NEXTSIMDIR/model/bin/nextsim.exec /usr/local/bin/nextsim.exec \
&&  ln -s $NEXTSIMDIR/model/run_in_docker.sh /usr/local/bin/run_in_docker.sh \
&&  echo $NEXTSIMDIR/lib/ >> /etc/ld.so.conf \
&&  ldconfig

WORKDIR /
ENTRYPOINT ["/usr/local/bin/run_in_docker.sh"]

