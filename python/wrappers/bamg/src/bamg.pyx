from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "bamg_wrap.hpp" namespace "PyWrap":
   vector[vector[double]] interpMeshToPointsCpp(
         vector[int] index,
         vector[double] xnods,
         vector[double] ynods,
         vector[vector[double]] data,
         vector[double] xout,
         vector[double] yout,
         bool isdefault, double defaultvalue)

# Import the Python-level symbols of numpy
# import numpy as np

# Import the C-level symbols of numpy
# cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
# np.import_array()


def interpMeshToPoints(index,xnods,ynods, data,xout,yout,
      isdefault=False, defaultvalue=1.e-24 ):

   if type(data)==type({}):
      # if data is a dictionary,
      #  convert to list for input
      interp_in   = []
      interp_out  = {}
      for key in data:
         interp_in.append(data[key])

      out   = interpMeshToPointsCpp(index,xnods,ynods, interp_in,xout,yout,
                  isdefault=False, defaultvalue=1.e-24 )

      # return output as dictionary
      for i,key in enumerate(data):
         interp_out.update({key:out[i]})

      return interp_out

   else:
      # if data is a list of vectors or 2d numpy array:
      out =  interpMeshToPointsCpp(index,xnods,ynods, data,xout,yout,
                  isdefault=False, defaultvalue=1.e-24 )

   if type(data)==type([]):
      # return list if list is input
      return list(out)
   else:
      return out
