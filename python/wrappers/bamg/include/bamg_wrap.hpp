#ifndef __BAMGWRAP_H
#define __BAMGWRAP_H 1

#include <vector>

namespace PyWrap
{

typedef std::vector<double> dbl_vec;
typedef std::vector<dbl_vec> dbl_vec2d;

dbl_vec2d interpMeshToPointsCpp(
      std::vector<int> index,
      dbl_vec xnods,
      dbl_vec ynods,
      std::vector< std::vector<double> > data,
      dbl_vec xout,
      dbl_vec yout,
      bool isdefault, double defaultvalue);

} //namespace PyWrappers

#endif
