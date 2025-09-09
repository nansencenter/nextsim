### Installation (that has worked with OSX 10.10) ###

#### Installation of the required compilers ####

1. Install MacPorts (`port`):
   * follow the instructions on the MacPorts [installation page](https://www.macports.org/install.php).
   * Check that `/opt/local/bin` and `/opt/local/sbin` are defined in PATH (`echo $PATH` in a new command window).

2. Install `gcc48` (or more recent versions) with MacPorts:
    ```
    sudo port search gcc
    sudo port install gcc48
    sudo port select --list gcc
    sudo port select --set gcc mp-gcc48
    ```

3. Install `openmpi-gcc48` (or a version compatible with your `gcc`) with MacPorts:
   ```
   sudo port search openmp
   sudo port install openmpi-gcc48
   sudo port select --list mpi
   sudo port select --set mpi openmpi-gcc48-fortran
   sudo ln -s /opt/local/bin/mpicxx /opt/local/bin/mpic++
   ```


#### Installation of the required libraries

1. Install the Boost C++ libraries
	1) Download version 1.55 of boost on http://www.boost.org (It is better to restart from here if you had an upgrade of your os).
	2) copy `bconfigure.sh` and `binstall.sh` from `/nextsim/scripts/boost/` to your boost directory.
	3) Add the line `using mpi ;` at the end of the file `tools/build/v2/user-config.jam` (This file could also be copied in your home directory for boost from the version 1.56).
	4) type the command `unset BOOST_DIR`
	5) check if `bconfigure.sh` corresponds to your architecture.
	6) From boost directory, type: `./bconfigure.sh`
	7) From boost directory, type: `./binstall.sh` (sudo password is required during the process).

2. Install NetCDF
	1) instal hdf5 via macport: `sudo port install hdf5`
	2) download latest stable C version of netcdf on [the NetCDF dowloads page](http://www.unidata.ucar.edu/downloads/netcdf/index.jsp) and unzip it.
	3) copy `configure_c.sh` from `/nextsim/scripts/netcdf`
	4) From netcdf directory, type: `./configure_c.sh` then `make` `make check` and finally `sudo make install`
	5) download latest stable C++ (`cxx`) version of netcdf on [the NetCDF dowloads page](http://www.unidata.ucar.edu/downloads/netcdf/index.jsp) and unzip it.
	6) copy `configure_cxx.sh` from `/nextsim/scripts/netcdf`
	7) From netcdf-cxx directory, type: `./configure_cxx.sh` then `make` `make check` and finally `sudo make install`

3. Install Gmsh from source
	1) Download the source code from [the Gmsh website](http://geuz.org/gmsh/).
	2) In the gmsh directory do
```
mkdir lib
cd lib
cmake -DDEFAULT=0 -DENABLE_BUILD_LIB=1 ..
make -j 32 lib
sudo make install/fast
```
4. Set the `PATH` environment variable correctly
	1) Add these lines to your `.bash_profile` for nextSIM in C++:
```
export NEXTSIMDIR=$HOME/Developer/nextsim/
export GMSH_DIR=/usr/local/

export NETCDF_DIR=/opt/local/netcdf-cxx

export BOOST_DIR=/opt/local/boost

export DYLD_LIBRARY_PATH="/opt/local/boost/lib"
export DYLD_LIBRARY_PATH=NEXTSIMDIR/lib:$DYLD_LIBRARY_PATH
```

#### Compile neXtSIM itself 

1. Open a new command window

2. Either modify the Makefile to have the right link to openmpi or do the following:
```
sudo rm -rf /opt/local/include/openmpi-mp
sudo ln -sf /opt/local/include/openmpi-gcc48 /opt/local/include/openmpi-mp
sudo rm -rf /opt/local/lib/openmpi-mp
sudo ln -sf /opt/local/lib/openmpi-gcc48 /opt/local/lib/openmpi-mp
```

3. Type `make` in nextsim
4. For the model application (run nextsim)
5. go to `nextsim/model`
6. Type `make`

### Run neXtSIM 

1. type `bin/nextsim.exec` from the `nextsim/model` directory.

###  Create documentation 
#### Run Doxygen with docker
```
cd [path to nextsim source]
docker run --rm -v $(pwd):/data -it hrektts/doxygen doxygen
```
This automatically uses the file Doxyfile (made originally with the doxygen gui)
  in the nextsim directory. Outputs go into a new folder `nextsim-doc`:
```
ls nextsim-doc/
html  latex  man  rtf  xml
```
You can then see the results with (eg.) `firefox nextsim-doc/html/index.html &`
#### Install Doxygen locally
* Ubuntu:
  `sudo apt-get install doxygen graphviz doxygen-gui`
* Mac OSX:
  * GUI: download .dmg image from docker website and open it to install 
  * command line: `sudo port install doxygen graphviz`
#### Run locally
* Ubuntu
  * GUI:
  `doxywizard &`
  * command line:
  ```
  cp Doxyfile Doxyfile-edited.cfg
  [edit new file and change entries for INPUT and OUTPUT_DIRECTORY]
  doxygen Doxyfile-edited.cfg
  ```
* Mac OSX:
  * GUI: open Doxygen.app file
  * command line: same as ubuntu
  * Note that Doxygen is not fully functional on Mac OSX (graphviz not working correctly) but GUI could still be used to create a config file
    to run with docker

