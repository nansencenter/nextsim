export BOOST_PREFIX=/opt/local/boost-1.60
#export BOOST_PREFIX=/opt/local/boost-1.66

./bootstrap.sh \
	--prefix=$BOOST_PREFIX \
	--without-libraries=python \
	--with-toolset=gcc
