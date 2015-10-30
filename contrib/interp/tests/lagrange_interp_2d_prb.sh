#!/bin/bash
#
g++ -c -I/$HOME/include lagrange_interp_2d_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling lagrange_interp_2d_prb.cpp"
  exit
fi
#
g++ lagrange_interp_2d_prb.o /$HOME/libcpp/$ARCH/lagrange_interp_2d.o \
                             /$HOME/libcpp/$ARCH/test_interp_2d.o \
                             /$HOME/libcpp/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lagrange_interp_2d_prb.o"
  exit
fi
#
rm lagrange_interp_2d_prb.o
#
mv a.out lagrange_interp_2d_prb
./lagrange_interp_2d_prb > lagrange_interp_2d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lagrange_interp_2d_prb."
  exit
fi
rm lagrange_interp_2d_prb
#
echo "Program output written to lagrange_interp_2d_prb_output.txt"
