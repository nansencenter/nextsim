FROM boost_petsc_gmsh:0.0.4

RUN apt-get update && apt-get install -y \
    libnetcdf-dev \
    libnetcdf-c++4-dev \
    ssh \
&& rm -rf /var/lib/apt/lists/*

ENV SRCDIR=/src

COPY Makefile $SRCDIR/
COPY contrib $SRCDIR/contrib
COPY core $SRCDIR/core
WORKDIR $SRCDIR
#RUN make

#COPY model $SRCDIR/model
#WORKDIR $SRCDIR/model
#RUN make

#COPY model/nextsim.sh /usr/local/bin/nextsim.sh && chmod a+x /usr/local/bin/nextsim.sh

WORKDIR $SRCDIR
#ENTRYPOINT ["/usr/local/bin/nextsim.sh"]
CMD [ "/bin/bash" ]
