export BOOST_PREFIX=/opt/local/boost-1.60
#export BOOST_PREFIX=/opt/local/boost-1.66

GCCDPC=$(echo `gcc -dumpversion | cut -f1-2 -d.` \>= 5.2 | sed -e 's/\.//g' | bc)
STDFLAG="\"-std=c++11\""
if [ "$GCCDPC" == "1" ]; then
	STDFLAG="\"-std=c++14\""
fi

sudo ./b2 -j8 -a\
	 --layout=tagged \
	 --debug-configuration \
	 --prefix=$BOOST_PREFIX \
	 toolset=gcc \
	 variant=release \
	 threading=single,multi \
	 link=shared,static \
	 cxxflags=$STDFLAG

sudo ./b2 install --prefix=$BOOST_PREFIX

# linkflags="-stdlib=libc++"
