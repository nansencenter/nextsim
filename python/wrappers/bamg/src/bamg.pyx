from libcpp.vector cimport vector
from libcpp cimport bool

# declare function to call
cdef extern from "bamg_wrap.hpp" namespace "PyWrap":
   vector[vector[double]] interpMeshToPointsCpp(
         vector[int] index,
         vector[double] xnods,
         vector[double] ynods,
         vector[vector[double]] data,
         vector[double] xout,
         vector[double] yout,
         bool isdefault, double defaultvalue)

def interpMeshToPoints(index,xnods,ynods, data,xout,yout,
      usedefault=False, defaultvalue=1.e-24 ):
   """
   call: bamg.interpMeshToPoints(index,xnods,ynods, data,xout,yout,
      isdefault=False, defaultvalue=1.e-24 )

   INPUTS:
   [numpy arrays] index,xnods,ynods: index (element->nodes map) and node coords

   data:
   - to be interpolated
   - possible types:
   1. list of numpy arrays, where each entry is a different variable
   2. 2d numpy array, where rows correspond to different variables
   3. dictionary of numpy arrays, where each key is a variable name

   [numpy arrays] xout,yout: coords of output location
   usedefault: True, then use defaultvalue input for points outside the mesh

   OUTPUTS:
   interpolated data
   - same type as 'data'

   Uses bamg function InterpMeshToMesh2dx

   """

   if type(data)==type({}):
      # if data is a dictionary,
      #  convert to list for input
      interp_in   = []
      interp_out  = {}
      for key in data:
         interp_in.append(data[key])

      out   = interpMeshToPointsCpp(index,xnods,ynods, interp_in,xout,yout,
                  usedefault, defaultvalue )

      # return output as dictionary
      for i,key in enumerate(data):
         interp_out.update({key:out[i]})

      return interp_out

   else:
      # if data is a list of vectors or 2d numpy array:
      out =  interpMeshToPointsCpp(index,xnods,ynods, data,xout,yout,
                  usedefault, defaultvalue )

   if type(data)==type([]):
      # return list if list is input
      return list(out)
   else:
      return out
