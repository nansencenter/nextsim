import os
import numpy as np
import nextsim_classes as nsc
import nextsim_funs as nsf
import nextsim_plot as nsp
import wim_grid_utils as wgu
plt   = nsp.plt

sddir = os.getenv('SIMDATADIR')
gmf   = sddir+'/data/mesh/wim_grid_FS_5km_split2.msh'

# get mesh object
gmsh  = nsc.gmsh_mesh(gmf)

# make a grid covering the same area
grid_params = gmsh.define_grid()
arrays      = wgu.make_grid(grid_params)
xout        = arrays['px'] # centre points of grid
yout        = arrays['py']

# define a land mask using nsf.interp_mesh_to_points
# (calls InterpMeshToMesh2dx)
data_in     = [np.zeros((gmsh.num_triangles))] #interp zeros to grid

# what to do outside the mesh
usedef   = True         # use the default value (defval) (if False, use the nearest point inside the mesh)
defval   = np.double(1) # any points outside the mesh get given the value 1.

# call the wrapper
land_mask   = nsf.interp_mesh_to_points(gmsh,xout,yout,data_in,
                  usedefault=usedef,defaultvalue=defval)[0]

plt.imshow(land_mask.transpose(),origin='lower')
plt.show()
