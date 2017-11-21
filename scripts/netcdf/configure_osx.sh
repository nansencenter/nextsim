export NETCDF_DIR=/opt/local/netcdf

LDFLAGS=-L/opt/local/lib \
./configure --prefix=$NETCDF_DIR \
	   --enable-netcdf-4 \
	   --enable-dynamic-loading \
	   --enable-extra-example-tests
