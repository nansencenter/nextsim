FROM akorosov/boost_petsc_gmsh:0.0.5

RUN apt-get update && apt-get install -yq --no-install-recommends \
    libnetcdf-dev \
    libnetcdf-c++4-dev \
    ssh \
    valgrind \
&&  apt-get clean \
&&  rm -rf /var/lib/apt/lists/*

ENV OPENMPI_INCLUDE_DIR=/usr/lib/x86_64-linux-gnu/openmpi/include \
    OPENMPI_LIB_DIR=/usr/lib/x86_64-linux-gnu/openmpi/lib \
    BOOST_INCDIR=/opt/local/boost/include \
    BOOST_LIBDIR=/opt/local/boost/lib \
    PETSC_DIR=/opt/local/petsc \
    GMSH_DIR=/opt/local/gmsh \
    PATH=$PATH:/src/model/bin \
    NEXTSIMDIR=/src \
    NEXTSIM_MESH_DIR=/mesh \
    NEXTSIM_DATA_DIR=/data


WORKDIR /src
