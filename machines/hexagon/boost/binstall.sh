export BOOST_PREFIX=$WORK/packages/boost

./b2 -j8 \
	--debug-configuration \
	--layout=tagged \
	--debug-configuration \
	--prefix=$BOOST_PREFIX \
	#toolset=darwin \
    #toolset=gcc \
	variant=release \
	threading=single,multi \
	link=shared,static \
	cxxflags="-std=c++11"

./b2 install --prefix=$BOOST_PREFIX

# linkflags="-stdlib=libc++"
