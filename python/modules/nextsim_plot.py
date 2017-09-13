from matplotlib import pyplot as plt


# ================================================================
class plot_object:
   def __init__(self):
      self.fig   = plt.figure()
      self.ax    = self.fig.add_subplot(1,1,1)
      return
# ================================================================


# ================================================================
def plot_mesh_data(mesh_obj,pobj=None,data=None,clabel=None,plot_grid=False,units='km',
      plot_coast=False,show=True,figname=None,colorbar=True):
   from matplotlib import patches,cm,collections
   import numpy as np

   if pobj is None:
      pobj  = plot_object()

   fig   = pobj.fig
   ax    = pobj.ax
   if units=='km':
      sfac  = 1e-3 #m to km
   elif units=='m':
      sfac  = 1
   else:
      raise ValueError('units should be km or m')

   # get nodes and indices
   nodes_x,nodes_y   = mesh_obj.get_nodes_xy()
   indices           = mesh_obj.get_indices("triangles",numbering='gmsh')
   Nn = mesh_obj.num_nodes
   Ne = mesh_obj.num_triangles
      
   # set axes ranges
   ax.set_xlim([sfac*nodes_x.min(),sfac*nodes_x.max()])
   ax.set_ylim([sfac*nodes_y.min(),sfac*nodes_y.max()])
   ax.set_aspect('equal')

   if data is None:
      # test data: zeros (to just plot mesh)
      zz       = np.zeros(Ne)
      colorbar = False
   elif len(data)!=Ne:
      # TODO scalars on nodes
      # TODO vectors on nodes - plot magnitude or component
      raise ValueError("plotting of nodal data not yet implemented")

   # get patches
   patch_list  = []
   for inds in indices:
      ccl   = []
      for n in inds:
         ccl.append((sfac*nodes_x[n],sfac*nodes_y[n]))
      ccl.append(ccl[0]) # close the contour
      patch_list.append(patches.Polygon(ccl,True))

   pc = collections.PatchCollection(patch_list, cmap=cm.jet, alpha=1)
   pc.set_array(data)

   if plot_grid:
      pc.set_edgecolor('k')
   else:
      pc.set_linewidth(0)

   ax.add_collection(pc)
   ax.set_xlabel('x, km')
   ax.set_ylabel('y, km')

   if colorbar:
      from mpl_toolkits.axes_grid1 import make_axes_locatable as MAL
      divider  = MAL(ax)
      cax   = divider.new_horizontal(pad=.45,size="5%",pack_start=False)
      fig.add_axes(cax)

      cbar  = fig.colorbar(pc,cax=cax,orientation="vertical")
      if clabel is not None:
         cbar.set_label(clabel,rotation=270)

   fig.tight_layout()

   
   if plot_coast and mesh_obj.boundary is not None:
      mesh_obj.boundary.plot(pobj=pobj,sort=True,show=False,units=units,thick_lines=plot_grid)

   if figname is not None:
      print("Saving as "+figname)
      fig.savefig(figname)

   if show:
      plt.show(fig)

   plt.close(fig)
   return pobj
# ================================================================
