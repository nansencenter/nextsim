sudo ./b2 -j8 \
	 --layout=tagged \
	 --debug-configuration \
	 --prefix=$BOOST_PREFIX \
	 toolset=darwin \
	 variant=release \
	 threading=single,multi \
	 link=shared,static \
	 cxxflags="-std=c++11"

sudo ./b2 install --prefix=$BOOST_PREFIX

# linkflags="-stdlib=libc++"
