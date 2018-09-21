# How to use Docker to build and run neXtSIM

## 1. Understand the concept

Docker is a platform for developers and sysadmins to develop, deploy, and run applications with containers.

Docker has two main concepts:
* image - a file with executables, binaries and environment
* container - a running instance of image

The workflow with Docker is to:
1. Build an *image* with the the nextsim exectuable and all needed libraries
2. Run a *container* with the application.

## 2. Install Docker

Firts, you need to install Docker Community Edition (Docker CE) on your platform.
See Docker overview and installation tips for Linux and Mac [here](https://docs.docker.com/install/)

Pay attention to [Linux optional post-installation steps](to https://docs.docker.com/install/linux/linux-postinstall/)

On macOS, once Docker is installed you need to go to the Docker icon and select Preferences->Advanced and set memory to at least 4 GB (the exact value needs testing). Otherwise the model won't compile.

## 3. Build an image

Image is built based on recipes in a `Dockerfile`. NeXtSIM repository contains Dockerfile for
compiling of bamx, mapx, core and the model code. This Dockerfile is based on another image
`boost_petsc_gmsh` available [here](https://github.com/nansencenter/docker-boost-petsc-gmsh)

#### 3.1 Build from scratch

To build an image you need to clone the repository, go to nextsim directory and run the following command:
```
docker build . -t nextsim
```
It will do the following:
* pull the Ubuntu image with Boost, PETSC and GMSH preinstalled
* spin up a temporary container based on this image
* copy the model source code into /opt/local/nextsim
* compile the code
* save the image with name `nextsim` (as provided with the `-t` parameter above)
* stop and remove the temporary container

#### 3.2 Recompile the code

Docker caches the images that are built. It means that if you run the above command again, without
changing the nextsim code, it will build the image from cache very fast.

If you want to introduce changes into the model, change the code and run `docker build . -t nextsim`
again. Docker will identify which files were changed and will run building of the image (including
compilation of the code) again.

## 4. Run the neXtSIM container

The image is built so that there are two options to run a container:
1. run `/bin/bash` and explore the environment inside container
2. run nextsim using MPI

#### 4.1 Run bash

If you only want to run /bin/bash, execute:
```
docker run --rm -it nextsim
```
The option `--rm` tells docker to remove the container after you exit from bash.
The option `-it` tells docker to run the container in foreground and provide interactive access to TTY.
`nextsim` is the name of the image that you have built previously.

#### 4.2 Run neXtSIM

If you want to run nextsim you also need to provide options to docker how to mount directories with
data, mesh and forecasts. By default, nextsim inside the container expects two directories:
```
export NEXTSIM_MESH_DIR=/mesh
export NEXTSIM_DATA_DIR=/data
```
Therefore you need to provide mounting of these two directories and a directory for output.

In addition you need to provide options to nextsim with the path to the config file and number of CPUs.
Note that all paths (to the config file, to the mesh, to the data, to outputs) are given in the
container file system.

An example command can look like the following:
```
docker run -it --rm \
    -v /home/user/nextsim/data/data_links:/data \
    -v /home/user/nextsim/mesh/mesh_links:/mesh \
    -v /home/user/output:/output \
    nextsim /output/test.cfg 7
```
It will, for example, mount directory `/home/user/nextsim/data/data_links` from the host as
`/data` on container.

The last two parameters (`/output/test.cfg` and `7`) are giving the location of the config file
and number of CPUs. Since we mount `/home/user/output` as `/output` the config file on a host
computer should be located in `/home/user/output/test.cfg`.

One more option `--security-opt seccomp=unconfined` is apparently needed to run MPI in container.

An example script to run model in a container can be found here:
[run_nextsim_container.sh](https://github.com/nansencenter/nextsim-env/blob/master/machines/tallinn/run_nextsim_container.sh)

#### 4.3 Debug neXtSIM

If you want to debug neXtSIM without rebuilind the entire image you can mount the nextsim
directory into /opt/local/nextsim inside the container
(with option `-v /path/on/the/host/nextsim:/opt/local/nextsim`). Then the compiled code will be
replaced with the source code, available for editing on your host machine. You should run the
container without specifying the config file and number of CPUs to enter bash inside the container.
Then you can compile the code (e.g. run `make All` from `/opt/local/nextsim`) and run the
model inside the container.

An example script to run container for debugging can be found here:
[run_nextsim_container_debug.sh](https://github.com/nansencenter/nextsim-env/blob/master/machines/tallinn/run_nextsim_container_debug.sh)
