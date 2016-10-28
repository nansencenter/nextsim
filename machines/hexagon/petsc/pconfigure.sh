export PETSC_PREFIX=$WORK/packages/petsc/3.6.3-mkl

./configure \
    -with-python \
    --with-debugging=0 \
    --with-c-support=1 \
    --with-c++-support=1 \
    --with-shared-libraries=1 \
    COPTFLAGS=-O3 CXXOPTFLAGS=-O3 FOPTFLAGS=-O3 \
    --with-mpi-dir=/opt/cray/mpt/7.2.5/gni/mpich2-gnu/51 \
    --with-blas-lapack-dir=$MKLROOT/lib/intel64 \
    --CFLAGS='-O3' \
    --CXXFLAGS='-O3' \
    --FFLAGS='-O3' \
    --with-blacs \
    --download-blacs=yes \
    --with-scalapack=1 \
    --download-scalapack=yes \
    --with-mumps=1 \
    --download-mumps=yes \
    --with-superlu=1 \
    --download-superlu=yes \
    --with-superlu_dist=1 \
    --download-superlu_dist=yes \
    --with-ml=1 \
    --download-ml=yes \
    --with-suitesparse=1 \
    --download-suitesparse=yes \
    --with-metis=1 \
    --download-metis=yes \
    --with-parmetis=1 \
    --download-parmetis=yes \
    --with-ptscotch=1 \
    --download-ptscotch=yes \
    --with-c2html=0 \
    --prefix=$PETSC_PREFIX \
    --with-mpiexec=aprun \

# --with-blas-lapack-lib="-L/opt/cray/libsci/13.2.0/GNU/4.9/x86_64/lib -lsci_gnu_mp" \
