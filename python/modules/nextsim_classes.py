# ================================================================================
class mesh_physical_name:
   def __init__(self,name,ident,topodim):
      self.name      = name
      self.ident     = ident
      self.topodim   = topodim
      return
# ================================================================================

# ================================================================================
class mesh_element:
   def __init__(self,ident,vertices,physical=None,elementary=None):
      self.ident        = ident
      self.vertices     = vertices
      self.num_vertices = len(vertices)
      self.physical     = physical
      self.elementary   = elementary
      return

   def get_coords(self,xnod,ynod):
      ccl  = []
      for n in self.vertices:
         ccl.append((xnod[n],ynod[n]))

      return ccl
# ================================================================================

# ================================================================================
class gmsh_boundary:

   # ==================================================================
   def __init__(self,exterior,islands=None,open_boundaries=None,coastal_boundaries=None):

      import numpy as np

      self.exterior_polygon   = exterior           # Shapely Polygon
      self.island_polygons    = islands            # List of shapely Polygons
      self.open_boundaries    = open_boundaries    # List of shapely LineStrings
      self.coastal_boundaries = coastal_boundaries # List of shapely LineStrings
      
      # get x-y range
      xe,ye       = exterior.exterior.xy
      self.xmin   = np.min(xe)
      self.xmax   = np.max(xe)
      self.ymin   = np.min(ye)
      self.ymax   = np.max(ye)

      # estimate resolution
      Ne    = len(xe)-1
      lens  = np.zeros((Ne))
      for i in range(Ne):
         x0 = xe[i]
         y0 = ye[i]
         x1 = xe[i+1]
         y1 = ye[i+1]
         lens[i]  = np.hypot(x1-x0,y1-y0)

      self.resolution   = np.mean(lens)

      return
   # ==================================================================


   # ==================================================================
   def iswet(self,x,y):
      """
      CALL: self.iswet(x,y)
      use matplotlib.path to test if multiple points are contained
      inside the polygon self.exterior_polygon
      INPUTS:
      x,y: numpy arrays with x and y coordinates to test respectively
      coords: list of coordinates of polygon's boundary (as tuples or lists)
      [(x0,y0),(x1,y1),...] or [[x0,y0],[x1,y1],...]
      RETURNS:
      mask of same shape as x and y
      > mask element is True/False if corresponding point is inside/outside the polygon
      """
      
      from matplotlib import path
      import numpy as np
      
      Nx    = x.size
      shp   = x.shape
      x     = x.reshape(Nx)
      y     = y.reshape(Nx)

      coords   = list(self.exterior_polygon.exterior.coords)
      bbPath   = path.Path(np.array(coords))

      # test if inside external polygon
      # points to test [(xi,yi)]
      # coords2  = [(x[i],y[i]) for i in range(Nx*Ny)]
      coords2  = np.array([x,y]).transpose()
      mask     = np.array(bbPath.contains_points(coords2),dtype=bool)

      for pp in self.island_polygons:
         # test not inside island polygons
         coords   = list(pp.exterior.coords)
         bbPath   = path.Path(np.array(coords))
         tmp      = np.array(bbPath.contains_points(coords2),dtype=bool)
         mask     = np.logical_and(mask,np.logical_not(tmp))

      return mask.reshape(shp)
   # ==================================================================


   # ==================================================================
   def define_grid(self,**kwargs):

      import nextsim_funs as nsf
      return nsf.define_grid(self,**kwargs)
   # ==================================================================


   # ==================================================================
   def plot(self,pobj=None,sort=True,show=True,units="km",thick_lines=False):

      import nextsim_plot as nsp
      import numpy as np
      plt   = nsp.plt

      if pobj is None:
         pobj  = nsp.plot_object()

      fig   = pobj.fig
      ax    = pobj.ax
      if units=="km":
         sfac  = 1e-3 #m to km
      elif units=="m":
         sfac  = 1
      else:
         raise ValueError('units should be km or m')

      if thick_lines:
         # lines need to stand out more sometimes
         lw = 5.35
      else:
         lw = 1.35

      # ===============================================================
      # plot exterior
      if not sort:
         pp    = self.exterior_polygon
         cc    = list(pp.exterior.coords)
         x,y   = np.array(cc).transpose()
         ax.plot(sfac*x,sfac*y,'k',linewidth=lw)
      else:
         if self.open_boundaries is not None:
            for ll in self.open_boundaries:
               cc    = list(ll.coords)
               x,y   = np.array(cc).transpose()
               ax.plot(sfac*x,sfac*y,'c',linewidth=lw)

         if self.coastal_boundaries is not None:
            for ll in self.coastal_boundaries:
               cc    = list(ll.coords)
               x,y   = np.array(cc).transpose()
               ax.plot(sfac*x,sfac*y,'k',linewidth=lw)
      # ===============================================================

      if self.island_polygons is not None:
         for pp in self.island_polygons:
            cc    = list(pp.exterior.coords)
            x,y   = np.array(cc).transpose()
            if not sort:
               ax.plot(sfac*x,sfac*y,'k',linewidth=lw)
            else:
               ax.plot(sfac*x,sfac*y,'r',linewidth=lw)

      if show:
         plt.show(fig)

      return pobj
   # ==================================================================

# ================================================================================

# ================================================================================
class gmsh_mesh:

   # ================================================================
   def __init__(self,meshfile,ordering='GMSH',mppfile=None,mapping=None):
      # TODO: add BAMG ordering permutation

      import nextsim_funs as nsf

      self.meshfile  = meshfile
      self.read_meshfile()

      # stereographic projection
      if mapping is not None:
         # projection already provided
         self.mppfile   = None
         self.mapping   = mapping
      else:
         # projection needs to be created
         self.mppfile   = mppfile
         self.mapping   = nsf.mppfile_to_pyproj(self.mppfile)

      self.stereographic_projection()

      self.get_boundary()

      self.xmin   = self.boundary.xmin
      self.ymin   = self.boundary.ymin
      self.xmax   = self.boundary.xmax
      self.ymax   = self.boundary.ymax

      self.resolution   = nsf.get_resolution(self)
      return
   # ================================================================

   # ================================================================
   def read_meshfile(self):

      import numpy as np
      import nextsim_funs as nsf

      self.version         = None
      self.format          = None
      self.size            = None
      self.PhysicalNames   = {}

      fid   = open(self.meshfile,'r')
      lin   = fid.readline()
      
      # ============================================================
      # see if GMSH mesh
      if "$MeshFormat" in lin:
         # =========================================================
         # mesh format and GMSH version
         lin                  = fid.readline()
         version,fmt,size     = lin.split()
         self.version         = float(version)
         self.format          = int(fmt)
         self.size            = int(size)
         lin                  = fid.readline()
         if "$EndMeshFormat" not in lin:
            raise ValueError("expecting $EndMeshFormat - not found")
         # =========================================================


         # =========================================================
         # GMSH meshes can have Physical Names defined at start
         self.num_physical_names = 0
         lin   = fid.readline()
         if "$PhysicalNames" in lin:
            lin   = fid.readline()
            self.num_physical_names = int(lin)

            # ========================================================================
            # loop over names
            for n in range(self.num_physical_names):
               lin                  = fid.readline()
               topodim,ident,name   = lin.split()
               if name[0]=='"':
                  name  = name[1:-1]
               self.PhysicalNames.update({int(ident):mesh_physical_name(name,int(ident),int(topodim))})
            # ========================================================================

            lin   = fid.readline()
            if "$EndPhysicalNames" not in lin:
               raise ValueError("expecting $PhysicalNames - not found")
            else:
               lin   = fid.readline()
         # =========================================================

      # =============================================================
      # get nodes
      nodstr   = lin.split()[0]
      if not ( ("$NOD" in lin) or ("$Nodes" in lin) or ("$ParametricNodes" in lin) ):
         raise ValueError("invalid nodes string '" + nodstr + "' in gmsh importer."\
               +" It should be either '$NOD','$Nodes' or '$ParametricNodes'.")

      lin               = fid.readline()
      self.num_nodes    = int(lin)
      self.nodes_id     = []
      self.nodes_map    = {}# map from node id to index
      self.nodes_lon    = []
      self.nodes_lat    = []
      for n in range(self.num_nodes):
         lin            = fid.readline()
         ident,c1,c2,c3 = lin.split()
         self.nodes_id.append(int(ident))
         self.nodes_map.update({int(ident):n})

         lon,lat  = nsf.xyz_to_lonlat(float(c1),float(c2),float(c3))
         self.nodes_lon.append(lon)
         self.nodes_lat.append(lat)

      lin   = fid.readline()
      if ("$End"+nodstr[1:] not in lin):
         raise ValueError("expecting $End"+nodstr[1:] +" - not found")

      lin   = fid.readline()
      self.nodes_lon  = np.array(self.nodes_lon)
      self.nodes_lat  = np.array(self.nodes_lat)
      # =============================================================


      if "$Elements" not in lin:
         raise ValueError("invalid elements string in gmsh importer.")

      self.element_num_vertices  = {1:2,2:3}

      lin               = fid.readline()
      self.num_elements = int(lin)
      self.edges        = []
      self.triangles    = []

      for n in range(self.num_elements):
         lin   = fid.readline()
         slin  = lin.split()

         elnumber       = int(slin[0])
         eltype         = int(slin[1])
         numtags        = int(slin[2])
         num_vertices   = self.element_num_vertices[eltype]

         jcpt  = 3
         for j in range(numtags):
            if j==0:
               physical = int(slin[jcpt])
            if j==1:
               elementary = int(slin[jcpt])
            jcpt  +=1

         # get vertices of elements
         verts = []
         for j in range(num_vertices):
            # convert from node id to node index
            # ind   = self.nodes_id.index(int(slin[jcpt]))
            ind   = self.nodes_map[int(slin[jcpt])]
            # ind   = int(slin[jcpt])-1 # faster, but less general
            verts.append(ind)
            jcpt  +=1

         if len(slin)!=jcpt:
            print(lin)
            raise ValueError("problem with number of vertices compared to element "+slin[0]+ " entry")

         elobj = mesh_element(elnumber,verts,physical=physical,elementary=elementary)
         if num_vertices==2:
            self.edges.append(elobj)
         elif num_vertices==3:
            self.triangles.append(elobj)

      lin   = fid.readline()
      if ("$EndElements" not in lin):
         raise ValueError("expecting $EndElements - not found")

      fid.close()

      self.num_edges       = len(self.edges)
      self.num_triangles   = len(self.triangles)
      # =============================================================

      return # read_meshfile
   # ================================================================


   # ================================================================
   def stereographic_projection(self):
      self.nodes_x,self.nodes_y  = self.mapping(self.nodes_lon,self.nodes_lat)
      return
   # ================================================================

   def get_nodes_xy(self):
      return self.nodes_x,self.nodes_y


   # ================================================================
   def get_indices(self,eltype="triangles",numbering='gmsh'):

      import numpy as np

      indices  = []
      if eltype=="triangles":
         els   = self.triangles
      else:
         els   = self.edges

      shift = 0
      if numbering=='bamg':
         # bamg nodes go from 1 to num_nodes
         # gmsh nodes go from 0 to num_nodes-1
         shift = 1

      for el in els:
         inds  = []
         for v in el.vertices:
            inds.append(v+shift)
         indices.append(inds)
         
      return np.array(indices)
   # ================================================================


   # ================================================================
   def get_boundary(self):
      import shapely.geometry as shg
      import shapely.ops as shops


      # ============================================================
      # get open and coast vals
      if len(self.PhysicalNames)==0:
         plabs = []
         cpt   = 0
         while (self.edges[cpt].physical not in plabs) and (len(plabs)<2):
            plabs.append(self.edges[cpt].physical)
            cpt  += 1
         self.coast_val = plabs[0]
         self.open_val  = plabs[1]
      else:
         for ident in self.PhysicalNames:
            if self.PhysicalNames[ident].name=="open":
               self.open_val = ident
            if self.PhysicalNames[ident].name=="coast":
               self.coast_val   = ident
      # ============================================================

      # ============================================================
      # sorting
      line_list   = [] # will be a list of shapely line strings
      open_list   = [] # will be a list of shapely line strings
      coast_list  = [] # will be a list of shapely line strings
      for edge in self.edges:
         ccl   = edge.get_coords(self.nodes_x,self.nodes_y)
         line_list.append(shg.LineString(ccl))
         if edge.physical==self.open_val:
            open_list.append(ccl)
         else:
            coast_list.append(ccl)
      # ===============================================================


      # ===============================================================
      # merge line segments to polygons
      lm    = shops.linemerge(line_list) #fewer multi-line strings
      if hasattr(lm,'geoms'):
         # multi line string
         polys = []
         for ll in lm:
            polys.append(shg.Polygon(ll))
      else:
         # single line string
         polys = [shg.Polygon(lm)]

      # sort by length to get exterior
      lengths  = []
      for pp in polys:
         lengths.append(len(pp.exterior.coords))

      if len(polys)==1:
         islands  = None
      else:
         lst      = sorted([(e,i) for i,e in enumerate(lengths)],reverse=True)
         polys    = [polys[i] for l,i in lst]
         islands  = polys[1:]

      exterior = polys[0]
      # ===============================================================


      # ===============================================================
      if len(open_list)==0:
         open_boundaries = None
      else:
         omm   = shops.linemerge(open_list) #fewer multi-line strings
         if hasattr(omm,'geoms'):
            open_boundaries = [e for e in omm]
         else:
            open_boundaries = [omm]

      # ===============================================================
      # here only give coasts that are not islands
      if len(coast_list)==0:
         coastal_boundaries = None
      else:
         cmm   = shops.linemerge(coast_list) #fewer multi-line strings
         coastal_boundaries = []
         if hasattr(cmm,'geoms'):
            for ll in cmm:
               if ll.intersects(exterior.exterior):
                  coastal_boundaries.append(ll)
         else:
            if cmm.intersects(exterior.exterior):
               coastal_boundaries = [cmm]

         if len(coastal_boundaries)==0:
            coastal_boundaries = None
      # ===============================================================

      self.boundary  = gmsh_boundary(exterior,islands=islands,\
                           open_boundaries=open_boundaries,coastal_boundaries=coastal_boundaries)

      return
   # ================================================================


   # ================================================================
   def plot_data(self,data=None,plotter="pyplot",**kwargs):

      indices  = self.get_indices("triangles")
      if data is None:
         import numpy as np
         data     = np.mean(self.nodes_lat[indices],axis=1)

      if plotter=="pyplot":
         # patches used in self.plot_data() is quicker than trying to
         # plot individual triangles
         import nextsim_plot as nsp
         clabel   = 'Latitude, degrees'
         kwargs.update({'plot_grid':True,'plot_coast':True})
         return nsp.plot_mesh_data(self,data=data,clabel=clabel,**kwargs)

      elif plotter=="mlab":
         # slower and can't zoom, but can rotate...
         from mayavi import mlab
         mlab.triangular_mesh(self.nodes_x,self.nodes_y,data,indices,\
               representation='wireframe',color=(0,0,0))
         mlab.view(0,0)
         mlab.show()

         return
   # ================================================================


   # ================================================================
   def plot_mesh(self,**kwargs):
      return self.plot_data(data=None,**kwargs)
   # ================================================================


   # ================================================================
   def plot_boundary(self,**kwargs):
      return self.boundary.plot(**kwargs)
   # ================================================================


# ===================================================================
class nextsim_mesh_info:
   def __init__(self,mesh_file,mapping=None,mppfile=None,gmsh_obj=None,gmsh_file=None,**kwargs):
      import nextsim_funs as nsf
      import numpy as np
      import os

      self.mesh_file       = os.path.abspath(mesh_file)
      self.mesh_file_info  = self.mesh_file.replace('.bin','.dat')

      self.variables,self.record_numbers,self.variable_types  = nsf.read_file_info(self.mesh_file_info)
      self.variable_lengths   = nsf.get_variable_lengths(self.mesh_file,self.variables,self.variable_types)

      self.num_nodes       = self.variable_lengths['Nodes_x']
      self.num_elements    = self.variable_lengths['Elements']/3
      self.num_triangles   = self.num_elements
      self.num_edges       = np.nan

      # stereographic projection
      if mapping is not None:
         # projection already provided
         self.mppfile   = None
         self.mapping   = mapping
      else:
         # projection needs to be created
         self.mppfile   = mppfile
         self.mapping   = nsf.mppfile_to_pyproj(self.mppfile)

      # boundary if given
      self.boundary  = None
      if gmsh_obj is not None:
         self.boundary  = gmsh_obj.boundary
      elif gmsh_file is not None:
         gmsh_obj = gmsh_mesh(gmsh_file,mapping=self.mapping,**kwargs)
         self.boundary  = gmsh_obj.boundary

      nodes_x,nodes_y   = self.get_nodes_xy()
      self.xmin         = np.min(nodes_x)
      self.ymin         = np.min(nodes_y)
      self.xmax         = np.max(nodes_x)
      self.ymax         = np.max(nodes_y)

      self.resolution   = nsf.get_resolution(self)
      return

   def get_var(self,vname):
      import nextsim_funs as nsf
      return nsf.get_array(vname,self.mesh_file,self.variables,self.variable_types)

   def get_vars(self,vname):
      import nextsim_funs as nsf
      return nsf.get_arrays(self.mesh_file,self.variables,self.variable_types)

   def get_nodes_xy(self):
      return self.get_var('Nodes_x'),self.get_var('Nodes_y')

   def get_indices(self,eltype="triangles",numbering='bamg',reshape=True):

      if eltype!="triangles":
         raise ValueError("nextsim_mesh_info object only has elements of type 'triangles'")

      shift = 0
      if numbering!='bamg':
         # bamg nodes go from 1 to num_nodes
         # gmsh nodes go from 0 to num_nodes-1
         shift = -1

      if not reshape:
         # return a vector [node_00,node_01,node_02,node_10,node_11,node_12,...]
         return self.get_var('Elements')+shift
      else:
         # reshape vector to 2d array out = [[node_00,node_01,node_02],[node_10,node_11,node_12],...],
         # where out[i] corresonds to i-th element
         return self.get_var('Elements').reshape((self.num_triangles,3))+shift

   def get_vertex_coords(self):
      vv       = self.get_vars()
      indices  = vv["Elements"].reshape((3,self.num_triangles)).transpose()
      verts    = []
      for inds in indices:
         x  = 1*vv["Nodes_x"][ind] # take a copy to destroy pointer
         y  = 1*vv["Nodes_y"][ind] # take a copy to destroy pointer
         verts.append((x,y))

   def plot(self,**kwargs):
      """
      plot the mesh
      self.plot(pobj=None,show=True,gmsh_boundary=None)
      """
      import nextsim_plot as nsp
      return nsp.plot_mesh_data(self,plot_grid=True,**kwargs)
# ================================================================================


# ================================================================================
class nextsim_binary_info:

   # =============================================================================
   def __init__(self,data_file,**kwargs):

      import os

      # full or relative path to file
      self.data_file       = os.path.abspath(data_file)
      self.data_file_info  = self.data_file.replace('.bin','.dat')
      ss             = self.data_file.split('/')
      self.basename  = ss[-1].strip('.bin')
      self.dirname   = '/' #self.data_file.strip('/'+ss[-1])
      for sdir in ss[:-1]:
         self.dirname  += (sdir+'/')

      # field info
      self._read_field_info()
      self.variables_all   = self.variables

      # mesh info
      ss          = self.dirname+'/'+self.basename.replace('field','mesh')
      mesh_file   = ss +'.bin'
      self.mesh_info = nextsim_mesh_info(mesh_file,**kwargs)

      #TODO read nextsim.log (mppfile to initialise projection, original mesh_file)
      #TODO or pass them in as arguments?

      return
   # =============================================================================

   # =============================================================================
   def _read_field_info(self):

      import nextsim_funs as nsf
      from datetime import datetime as DT,timedelta as TD

      # get the variables and their record no
      self.variables,self.record_numbers,self.variable_types  = nsf.read_file_info(self.data_file_info)

      Time  =  nsf.get_array("Time",self.data_file,self.variables,self.variable_types)[0]
      
      refdate        = DT(1900,1,1)
      self.datetime  = refdate+TD(Time)
      self.datetimes = [self.datetime]

      return
   # =============================================================================

   # =============================================================================
   def get_var(self,vname):
      import nextsim_funs as nsf
      return nsf.get_array(vname,self.data_file,self.variables,self.variable_types)

   def get_vars(self):
      import nextsim_funs as nsf
      return nsf.get_arrays(self.data_file,self.variables,self.variable_types)

   def plot_mesh(self,**kwargs):
      return self.mesh_info.plot(**kwargs)

   def plot_var(self,vname,**kwargs):
      import nextsim_plot as nsp
      return nsp.plot_mesh_data(self.mesh_info,data=self.get_var(vname),**kwargs)
   # =============================================================================

# ================================================================================


