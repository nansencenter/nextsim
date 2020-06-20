#!/bin/bash
#
cp pwl_interp_2d_scattered.hpp /$HOME/include
#
g++ -c -I /$HOME/include pwl_interp_2d_scattered.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling pwl_interp_2d_scattered.cpp"
  exit
fi
#
mv pwl_interp_2d_scattered.o ~/libcpp/$ARCH/pwl_interp_2d_scattered.o
#
echo "Library installed as ~/libcpp/$ARCH/pwl_interp_2d_scattered.o"
