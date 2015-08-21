/*This routine only used by Intel compler: */

#include "isnan.h"

#if defined(_HAVE_ADOLC_) && !defined(_WRAPPERS_)
template <> int xIsNan<adouble> (const adouble& X){
  return isnan(X.getValue());
}
#endif
