FROM akorosov/boost_petsc_gmsh:0.0.5

RUN apt-get update && apt-get install -y \
    libnetcdf-dev \
    libnetcdf-c++4-dev \
    ssh \
&& rm -rf /var/lib/apt/lists/*

ENV NEXTSIMDIR=/opt/local/nextsim
SHELL ["/bin/bash", "-c"]

COPY contrib/mapx $NEXTSIMDIR/contrib/mapx
WORKDIR $NEXTSIMDIR/contrib/mapx/src
RUN source /root/.nextsimrc && make

COPY contrib/bamg $NEXTSIMDIR/contrib/bamg
WORKDIR $NEXTSIMDIR/contrib/bamg/src
RUN source /root/.nextsimrc && make

COPY core $NEXTSIMDIR/core
WORKDIR $NEXTSIMDIR/core/src
RUN source /root/.nextsimrc && make

COPY model $NEXTSIMDIR/model
WORKDIR $NEXTSIMDIR/model
RUN source /root/.nextsimrc && make

RUN ln -s $NEXTSIMDIR/model/bin/nextsim.exec /usr/local/bin/nextsim.exec \
&&  ln -s $NEXTSIMDIR/model/run_in_docker.sh /usr/local/bin/run_in_docker.sh \
&&  echo $NEXTSIMDIR/lib/ >> /etc/ld.so.conf \
&&  ldconfig

WORKDIR /
ENTRYPOINT ["/usr/local/bin/run_in_docker.sh"]
#CMD [ "/bin/bash" ]
