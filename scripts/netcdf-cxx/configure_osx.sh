export NETCDF_DIR=/opt/local/netcdf
export NETCDFCXX_DIR=/opt/local/netcdf-cxx

export CPPFLAGS=-I$NETCDF_DIR/include
export LDFLAGS=-L$NETCDF_DIR/lib

./configure --prefix=$NETCDFCXX_DIR
