# mapx
cd contrib/mapx/src
make mrproper
make -j 4
make clean

# oasis coupling module
cd ../../../modules/oasis/src
make mrproper
make
make clean

# bamg
cd ../../../contrib/bamg/src
make mrproper
make -j 4
make clean

# gmsh 
cd ../../../core/src
make mrproper
make -j 4
make clean

# build nextsim.exec
cd ../../model/
make mrproper
make
make clean

