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

On macOS, once Docker is installed you need to go to the Docker icon and select
Preferences->Advanced and set memory to at least 4 GB (the exact value needs testing).
Otherwise the model won't compile.

## 3. Build an image

The Docker image is built based on recipes in a `Dockerfile`. NeXtSIM repository contains a
Dockerfile for installation of compilers and libraries. This Dockerfile is based on another image
`boost_petsc_gmsh` available [here](https://github.com/nansencenter/docker-boost-petsc-gmsh)

To build an image you need to clone the repository, go to nextsim directory and run the following
command:
```
docker build . -t nextsim
```
It will pull the Docker image with Ubuntu, Boost, PETSC and GMSH and install gcc, libnetcdf and
compile BAMG and MAPX libraries.

## 4. Compile the code

Now you can use this image for compiling nextsim core and model. You need to run a container:
```
docker run --rm -v `pwd`:/nextsim nextsim make docker -j8
```
These commands will:
* start a container whith all the required libraries
* set environment variables in the container (PATH, NEXTSIMDIR, etc..)
* mount the current folder (with nextsim source code) into /nextsim in container
* compile the code in the core and model directories and save the objects and executable in the
mounted directories. NB: the generated binary files will be available both for the host
(in the current directory) and for the container (in /nextsim).

If you want to recompile only the model code :
```
docker run --rm -v `pwd`:/nextsim nextsim make cleandocker
docker run --rm -v `pwd`:/nextsim nextsim make docker -j8
```

## 4. Run the neXtSIM executable inside a container

The image is built to run any executable inside a container. For example, if you only want to run
bash:
```
docker run -it --rm -v `pwd`:/nextsim nextsim bash
```
The option `--rm` tells docker to remove the container after you exit from bash.
The option `-it` tells docker to run the container in foreground and provide interactive access to TTY.
`nextsim` is the name of the image that you have built previously. `bash` - is the command ro run.

If you want to run nextsim you also need to provide options to docker how to mount directories with
data, mesh and outputs. Nextsim inside the container expects two directories:
```
NEXTSIM_MESH_DIR=/mesh
NEXTSIM_DATA_DIR=/data

```
Therefore you need to provide mounting of these two directories and a directory for output.
In addition you need to provide the command to run mpirun with the path to the config file, number
of CPUs and so on. Note that all paths (to the config file, to the mesh, to the data, to outputs)
are given in the container file system.

An example command can look like the following:
```
docker run -it --rm \
    --security-opt seccomp=unconfined \
    -v /home/user/nextsim:/nextsim
    -v /Data/sim/data:/data \
    -v /home/user/nextsim/mesh:/mesh \
    -v /home/user/output:/output \
    nextsim \
    mpirun --allow-run-as-root -np 8 nextsim.exec -mat_mumps_icntl_23 200 --config-files=/output/test.cfg
```
This example will mount the following directories:
* `/home/user/nextsim` on the host with source code and compiled binaries as `/nextsim` in container
* `/Data/sim/data` with all input data as `/data` in container
* `/home/user/nextsim/mesh` with mpp files and links to meshes as `/mesh` in container
* `/home/user/output` with model output as `/output` in container
If your directories (e.g. `/home/user/nextsim/mesh`) contain not the files but symbolic links to files,
you also need to mount the directories where the files are actually residing
(e.g. `-v /Data/sim/data:/Data/sim/data`)

One more option `--security-opt seccomp=unconfined` is apparently needed to run MPI in container.

An example script to run model in a container can be found here:
[run_nextsim_container.sh](https://github.com/nansencenter/nextsim-env/blob/master/machines/maud_antonk/run_nextsim_container.sh)

