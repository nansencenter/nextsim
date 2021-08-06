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
Dockerfile for compiling the model code. This Dockerfile is based on another image
`nextsim_base` openly available at Docker hub. The Dockerfile for nextsim_base is located in
[nextsim-env](https://github.com/nansencenter/nextsim-env) repo.

To build an image you need to clone the repository, go to nextsim directory and run the following
command:
```
docker build . -t nextsim
```
It will pull the Docker image with Ubuntu, Boost and GMSH and
compile BAMG and MAPX libraries and compile the core and the model itself.

Now you can run the model code which is inside the image:
```
docker run --rm nextsim nextsim.exec
```
This will, of course, crash because no configfile is given to the model. But it
shows that the model is already compiled and stored inside the image. In the section 4 below see
how to run the model correctly.

## 4. Run the neXtSIM executable inside a container

The image is built to run any executable inside a container. For example, if you only want to run
bash:
```
docker run -it --rm nextsim bash
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
These directories should contain either files or links to files. There are two helper scripts that
can create these links automatically. You can run these scripts inside the `nextsim` image:
```
# link meshes
docker run --rm \
    -v /Data/sim/:/Data/sim/ \
    -v /home/user/nextsim/mesh:/mesh \
    nextsim \
    link_meshes.sh /Data/sim/data/mesh /mesh

# link datafiles
docker run --rm \
    -v /Data/sim/:/Data/sim/ \
    -v /Data/nextsimf/:/Data/nextsimf/ \
    -v /home/user/nextsim/data:/data \
    nextsim \
    link_data.sh /Data/sim/data /Data/nextsimf/data /data
```

In addition you need to provide the command to run mpirun with the path to the config file, number
of CPUs and so on. Note that all paths (to the config file, to the mesh, to the data, to outputs)
are given in the container file system.

An example command can look like the following:
```
docker run -it --rm \
    --security-opt seccomp=unconfined \
    -v /home/user/nextsim/data:/data \
    -v /home/user/nextsim/mesh:/mesh \
    -v /home/user/output:/output \
    nextsim \
    mpirun --allow-run-as-root -np 8 nextsim.exec -mat_mumps_icntl_23 200 --config-files=/output/test.cfg
```
This example will mount the following directories:
* `/home/user/nextsim/data` with links to all input data as `/data` in container
* `/home/user/nextsim/mesh` with mpp files and links to meshes as `/mesh` in container
* `/home/user/output` with model output as `/output` in container
If your directories (e.g. `/home/user/nextsim/mesh`) contain not the files but symbolic links to files,
you also need to mount the directories where the files are actually residing
(e.g. `-v /Data/sim/data:/Data/sim/data`)

One more option `--security-opt seccomp=unconfined` is apparently needed to run MPI in container.

An example script to run model in a container can be found here:
[run_nextsim_container.sh](https://github.com/nansencenter/nextsim-env/blob/master/machines/maud_antonk/run_nextsim.sh)

## 5. It is still possible to compile the model code in-place. "In-place" means
that the binary objects and exceutbale are created in the directory on the host machine. For
compiling inplace you need to mount the current folder into `/nextsim` and set the working
directory with `-w` options:

```
docker run --rm -it -v /home/user/nextsim:/nextsim -w /nextsim/model nextsim make
```

**Remember:** if you want to run the code compiled in-place, mount the `/nextsim` directory.
