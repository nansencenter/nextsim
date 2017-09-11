#ifndef __BAMGWRAP_H
#define __BAMGWRAP_H 1

#include <vector>

namespace PyWrap
{

std::vector<std::vector<double>> interpMeshToPointsCpp(
      std::vector<int> index,
      std::vector<double> xnods,
      std::vector<double> ynods,
      std::vector<std::vector<double>> data,
      std::vector<double> xout,
      std::vector<double> yout,
      bool isdefault, double defaultvalue);

} //namespace PyWrappers

#endif
