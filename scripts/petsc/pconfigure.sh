export PETSC_PREFIX=/opt/local/petsc

./configure \
	-with-python \
	--with-cxx=mpicxx \
	--with-cc=mpicc \
	--with-fc=mpif90 \
	--with-mpiexec=mpiexec \
	--with-debugging=0 \
	--with-c-support=1 \
	--with-c++-support=1 \
	--with-shared-libraries=1 \
	--with-blacs \
	--download-blacs=yes \
	--with-parmetis=1 \
	--download-parmetis=yes \
	--with-scalapack=1 \
	--download-scalapack=yes \
	--with-mumps=1 \
	--download-mumps=yes \
	--with-suitesparse=1 \
	--download-suitesparse=yes \
	--with-c2html=0 \
	--with-metis=1 \
	--download-metis=yes \
	--prefix=$PETSC_PREFIX
