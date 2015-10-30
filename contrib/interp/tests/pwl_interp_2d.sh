#!/bin/bash
#
cp pwl_interp_2d.hpp /$HOME/include
#
g++ -c -I/$HOME/include pwl_interp_2d.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling pwl_interp_2d.cpp"
  exit
fi
#
mv pwl_interp_2d.o ~/libcpp/$ARCH/pwl_interp_2d.o
#
echo "Library installed as ~/libcpp/$ARCH/pwl_interp_2d.o"
