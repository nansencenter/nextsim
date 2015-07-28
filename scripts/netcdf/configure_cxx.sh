export NETCDF_DIR=/opt/local/netcdf
export NETCDFCXX_DIR=/opt/local/netcdf-cxx

CPPFLAGS=-I$NETCDF_DIR/include LDFLAGS=-L$NETCDF_DIR/lib \
		./configure --prefix=$NETCDFCXX_DIR
