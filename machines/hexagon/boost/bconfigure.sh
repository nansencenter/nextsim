export BOOST_PREFIX=$WORK/packages/boost

./bootstrap.sh \
	--prefix=$BOOST_PREFIX \
	--without-libraries=python \
	cxxflags="-arch i386 -arch x86_64" \
	address-model=32_64 \
	threading=single,multi
	#macos-version=10.10
