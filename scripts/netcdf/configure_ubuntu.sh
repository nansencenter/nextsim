export NETCDF_DIR=/opt/local/netcdf

export LDFLAGS="-L/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu/hdf5/serial"

export CPPFLAGS="-I/usr/lib/x86_64-linux-gnu/hdf5/serial/include"
./configure --prefix=$NETCDF_DIR \
	   --enable-netcdf-4 \
	   --enable-dynamic-loading \
	   --enable-extra-example-tests
