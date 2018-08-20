FROM boost_petsc_gmsh:0.0.4

RUN apt-get update && apt-get install -y \
    libnetcdf-dev \
    libnetcdf-c++4-dev \
    ssh \
&& rm -rf /var/lib/apt/lists/*

ENV NEXTSIMDIR=/src
SHELL ["/bin/bash", "-c"]

COPY Makefile $NEXTSIMDIR/
COPY contrib $NEXTSIMDIR/contrib
COPY core $NEXTSIMDIR/core
WORKDIR $NEXTSIMDIR/
RUN source /root/.nextsimrc && make

COPY model $SRCDIR/model
WORKDIR $SRCDIR/model
RUN source /root/.nextsimrc && make

#COPY model/nextsim.sh /usr/local/bin/nextsim.sh
#RUN mkdir -p /opt/local/nextsim/lib && mv $SRCDIR/lib/* /opt/local/nextsim/lib/
#RUN mv $SRCDIR/model/bin/nextsim.exec /usr/local/bin/ && chmod a+x /usr/local/bin/nextsim.sh
#RUN echo '/opt/local/nextsim/lib/' >> /etc/ld.so.conf \
#&&  ldconfig \

WORKDIR $NEXTSIMDIR
#ENTRYPOINT ["/usr/local/bin/nextsim.sh"]
CMD [ "/bin/bash" ]
