# find closest point on WIM grid,
# corresponding to a given (x,y) or (lon,lat) coordinate
import fns_get_data as Fdat
import numpy as np

# res   = Fdat.wim_results('out_cpp')
res   = Fdat.wim_results('.')
G     = res.get_grid()
X     = G['X']
Y     = G['Y']

if 1:
   # give x,y
   x,y   = 1064294.30630935,-920004.063012927 

R  = np.hypot(X-x,Y-y)
I  = np.where(R==R.min())
print(I)
print(x,X[I])
print(y,Y[I])
