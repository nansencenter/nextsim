export PETSC_PREFIX=/opt/local/petsc/3.6.3

./configure \
    -with-python \
    --with-debugging=0 \
    --with-c-support=1 \
    --with-c++-support=1 \
    --with-shared-libraries=1 \
    --with-cxx=mpicxx \
    --with-cc=mpicc \
    --with-fc=mpif90 \
    --CFLAGS='-O3' \
    --CXXFLAGS='-O3' \
    --FFLAGS='-O3' \
    --with-mpiexec=mpiexec \
    --with-blacs \
    --download-blacs=yes \
    --with-scalapack=1 \
    --download-scalapack=yes \
    --with-mumps=1 \
    --download-mumps=yes \
    --with-pastix=1 \
    --download-pastix=yes \
    --with-superlu=1 \
    --download-superlu=yes \
    --with-suitesparse=1 \
    --download-suitesparse=yes \
    --with-metis=1 \
    --download-metis=yes \
    --with-parmetis=1 \
    --download-parmetis=yes \
    --with-ptscotch=1 \
    --download-ptscotch=yes \
    --with-c2html=0 \
    --prefix=$PETSC_PREFIX



# --COPTFLAGS='-O3' \
# --CXXOPTFLAGS='-O3' \
# --FOPTFLAGS='-O3' \
# --FFLAGS='-mavx -mfma4 -O3 -march=native -mtune=native -Wa,-q' \
# --FOPTFLAGS='-mavx -mfma4 -O3 -march=native -mtune=native -Wa,-q' \
