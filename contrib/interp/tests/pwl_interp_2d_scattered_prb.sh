#!/bin/bash
#
g++ -c -I/$HOME/include pwl_interp_2d_scattered_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling pwl_interp_2d_scattered_prb.cpp"
  exit
fi
#
g++ pwl_interp_2d_scattered_prb.o /$HOME/libcpp/$ARCH/pwl_interp_2d_scattered.o \
                                  /$HOME/libcpp/$ARCH/test_interp_2d.o \
                                  /$HOME/libcpp/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pwl_interp_2d_scattered_prb.o."
  exit
fi
#
rm pwl_interp_2d_scattered_prb.o
#
mv a.out pwl_interp_2d_scattered_prb
./pwl_interp_2d_scattered_prb > pwl_interp_2d_scattered_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pwl_interp_2d_scattered_prb."
  exit
fi
rm pwl_interp_2d_scattered_prb
#
echo "Program output written to pwl_interp_2d_scattered_prb_output.txt"
