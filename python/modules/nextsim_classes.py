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

      print('Reading '+meshfile+'...\n')
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

      print('Finished reading '+meshfile+'.\n')
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

   def get_nodes_xy(self,dtype='float'):
      import numpy as np
      return np.array(self.nodes_x,dtype=dtype),\
               np.array(self.nodes_y,dtype=dtype)


   # ================================================================
   def get_indices(self,eltype="triangles",numbering='gmsh',asvector=True):

      import numpy as np

      indices  = []
      if eltype=="triangles":
         els      = self.triangles
         nverts   = 3
      else:
         els      = self.edges
         nverts   = 2

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

      if asvector:
         return np.array(indices).reshape((len(els)*nverts))
      else:
         # 
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

      indices  = self.get_indices(eltype="triangles",
            numbering='gmsh',asvector=False)

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


   # ================================================================
   def define_grid(self,**kwargs):
      import nextsim_funs as nsf
      return nsf.define_grid(self,**kwargs)
   # ================================================================


# ===================================================================
class nextsim_mesh_info:

   # ================================================================
   def __init__(self,mesh_file,mapping=None,mppfile=None,\
         boundary=None,gmsh_obj=None,gmsh_file=None,**kwargs):
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
      self.boundary  = boundary
      if boundary is None:
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
   # ================================================================


   # ================================================================
   def get_var(self,vname):
      import nextsim_funs as nsf
      return nsf.get_array(vname,self.mesh_file,self.variables,self.variable_types)
   # ================================================================


   # ================================================================
   def get_vars(self,vname):
      import nextsim_funs as nsf
      return nsf.get_arrays(self.mesh_file,self.variables,self.variable_types)
   # ================================================================


   # ================================================================
   def get_nodes_xy(self,dtype='float'):
      import numpy as np
      return np.array(self.get_var('Nodes_x'),dtype=dtype),\
               np.array(self.get_var('Nodes_y'),dtype=dtype)
   # ================================================================


   # ================================================================
   def get_indices(self,eltype="triangles",numbering='bamg',asvector=False):

      if eltype!="triangles":
         raise ValueError("nextsim_mesh_info object only has elements of type 'triangles'")

      shift = 0
      if numbering!='bamg':
         # bamg nodes go from 1 to num_nodes
         # gmsh nodes go from 0 to num_nodes-1
         shift = -1

      if asvector:
         # return a vector [node_00,node_01,node_02,node_10,node_11,node_12,...]
         return self.get_var('Elements')+shift
      else:
         # reshape vector to 2d array out = [[node_00,node_01,node_02],[node_10,node_11,node_12],...],
         # where out[i] corresonds to i-th element
         return self.get_var('Elements').reshape((self.num_triangles,3))+shift
   # ================================================================


   # ================================================================
   def get_vertex_coords(self):
      vv       = self.get_vars()
      indices  = vv["Elements"].reshape((self.num_triangles,3))
      verts    = []

      for inds in indices:
         x  = 1*vv["Nodes_x"][ind] # take a copy to destroy pointer
         y  = 1*vv["Nodes_y"][ind] # take a copy to destroy pointer
         verts.append((x,y))

      return verts
   # ================================================================


   # ================================================================
   def plot(self,**kwargs):
      """
      plot the mesh
      self.plot(pobj=None,show=True,gmsh_boundary=None)
      """
      import nextsim_plot as nsp
      return nsp.plot_mesh_data(self,plot_grid=True,**kwargs)
   # ================================================================


   # ================================================================
   def define_grid(self,**kwargs):
      import nextsim_funs as nsf
      return nsf.define_grid(self,**kwargs)
   # ================================================================

# gmsh_mesh class
# ================================================================================


# ================================================================================
class nextsim_binary_info:

   # =============================================================================
   def __init__(self,data_file,options=None,logfile=None,**kwargs):
      """
      nbi   = nextsim_classes.nextsim_binary_info(data_file,logfile=None,**kwargs)
      * data_file: eg field_0.bin
      * logfile: eg nextsim.log
      * kwargs are passed to nextsim_classes.nextsim_mesh_info()
      """

      import os

      # =======================================================
      # full or relative path to file
      self.data_file       = os.path.abspath(data_file)
      self.data_file_info  = self.data_file.replace('.bin','.dat')
      ss                   = self.data_file.split('/')
      self.basename        = ss[-1].strip('.bin')
      self.dirname         = '/' #self.data_file.strip('/'+ss[-1])
      for sdir in ss[:-1]:
         self.dirname  += (sdir+'/')
      # =======================================================


      # =======================================================
      if options is not None:
         self.logfile   = logfile
         self.options   = options
      else:
         # read log file to get options
         # - update kwargs if possible
         self.logfile   = logfile
         self._read_log()

      self._check_options(kwargs)
      # =======================================================


      # =======================================================
      # field info
      self._read_field_info()
      self.variables_all   = self.variables

      # mesh info
      ss          = self.dirname+'/'+self.basename.replace('field','mesh')
      mesh_file   = ss +'.bin'
      self.mesh_info = nextsim_mesh_info(mesh_file,**kwargs)
      # =======================================================

      return
   # =============================================================================


   # =============================================================================
   def _read_log(self):
      # reads log if present
      import nextsim_funs as nsf
      import os

      # =======================================================
      # read log file to get options
      if self.logfile is None:
         lf = self.dirname+'/nextsim.log'
         if os.path.exists(lf):
            self.logfile   = lf

      self.options   = nsf.read_nextsim_log(self.logfile)
      return
   # =============================================================================


   # =============================================================================
   def _check_options(self,kwargs):
      # checks options and updates kwargs to pass to nextsim_mesh_info

      if self.options is not None:
         import os

         # =======================================================
         # determine original mesh file if needed
         haveinfo = False
         for key in ['gmsh_file','gmsh_obj']:
            if key in kwargs:
               if key is not None:
                  haveinfo = True

         if not haveinfo:
            # can give extra info to the mesh object about the original mesh
            # (this gives the coastlines)
            gmsh_dirs   = {'simdatadir':os.getenv('SIMDATADIR')+'/data/mesh',\
                           'nextsimdir':os.getenv('NEXTSIMDIR')+'/mesh'}
            gmsh_file   = self.options['program_options']['simul']['mesh_filename']
            gmsh_dir    = self.options['program_options']['simul']['mesh_path']
            if gmsh_dir in gmsh_dirs:
               gmsh_dir = gmsh_dirs[gmsh_dir]
            kwargs.update({'gmsh_file':gmsh_dir+'/'+gmsh_file})
         # =======================================================

         # =======================================================
         # determine mapping if needed
         haveinfo = False
         for key in ['mppfile','mapping']:
            if key in kwargs:
               if key is not None:
                  haveinfo = True

         if not haveinfo:
            # can give extra info to the mesh object about the original mesh
            # (this gives the coastlines)
            mpp_dir  = os.getenv('NEXTSIMDIR')+'/data'
            mpp_file = mpp_dir+'/'+self.options['program_options']['simul']['proj_filename']
            kwargs.update({'mppfile':mpp_file})
         # =======================================================

      return
      # =======================================================


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
   # =============================================================================


   # =============================================================================
   def get_vars(self,vlist=None):
      import nextsim_funs as nsf
      if vlist is None:
         vlist = self.variables
      return nsf.get_arrays(self.data_file,vlist,self.variable_types)
   # =============================================================================


   # =============================================================================
   def plot_mesh(self,**kwargs):
      return self.mesh_info.plot(**kwargs)
   # =============================================================================


   # =============================================================================
   def plot_var(self,vname,**kwargs):
      import nextsim_plot as nsp
      return nsp.plot_mesh_data(self.mesh_info,data=self.get_var(vname),**kwargs)
   # =============================================================================


   # =============================================================================
   def interp_var(self,vlist,xout,yout,**kwargs):
      """
      self.interp_var(vlist,x_out,y_out,**kwargs)
      
      Inputs:

      vlist: list of variables
      x_out,y_out: lists/vectors with x,y coordinates
      kwargs (for nextsim_funs.interp_mesh_to_points):
         use_default=True - use default_value outside the mesh; if False, use the nearest point inside the mesh
         default_value=0., see use_default=True
      
      Output:
      dictionary with requested fields interpolated onto x_out,y_out
      """

      import nextsim_funs as nsf
      data_in  = {}
      for vname in vlist:
         data_in.update({vname:self.get_var(vname)})

      # do interpolation
      return nsf.interp_mesh_to_points(self.mesh_info,xout,yout,data_in,**kwargs)
   # =============================================================================


   # =============================================================================
   def imshow(self,vname):

      import numpy as np
      from matplotlib import pyplot as plt
      import wim_grid_utils as wgu

      # make a grid covering the same area as the mesh
      grid_params = self.mesh_info.define_grid()
      arrays      = wgu.make_grid(grid_params)
      xout        = arrays['px'] # centre points of grid
      yout        = arrays['py']

      # interpolate vname
      usedef   = True
      defval   = np.NaN
      zout     = self.interp_var([vname],xout,yout,use_default=usedef,default_value=defval)

      fig   = plt.figure()
      ax    = fig.add_subplot(111)
      im    = ax.imshow(zout[vname].transpose(),origin='upper')
      fig.colorbar(im)
      plt.show(fig)

      return
   # =============================================================================

# ================================================================================


# ================================================================================
class file_list:
   def __init__(self,rootdir='.',name_filter=None,\
         step1=None,step2=None,date1=None,date2=None,**kwargs):

      import os
      import numpy as np
      import nextsim_funs as nsf
      import time
      t0 = time.time()

      lst   = os.listdir(rootdir)

      # check the files
      filelist  = []
      steplist  = []
      datelist  = []
      for f in lst:

         # =============================================
         # check if we should include the file
         prefix   = 'field_'
         if name_filter is not None:
            prefix  += name_filter+'_'

         if prefix in f and '.bin' in f:
            # check if step is a number
            cstep = f[len(prefix):len(f)-4]
            try:
               step  = int(cstep)
            except:
               continue

            # check if it is in the step range
            if step1 is not None:
               date1 = None
               if step<step1:
                  continue

            if step2 is not None:
               date2 = None
               if step>step2:
                  continue

            ff    = os.path.abspath(rootdir)+'/'+f
            fdat  = ff.replace('.bin','.dat')
            dto   = nsf.get_nextsim_time(fdat)

            if type(date1) != type(None):
               if dto<date1:
                  continue

            if type(date2) != type(None):
               if dto>date2:
                  continue

            steplist.append(step)
            filelist.append(ff)
            datelist.append(dto)
         # finished sorting files
         # =============================================

      # ======================================================================
      # now sort the times into order
      lst   = sorted([(e,i) for i,e in enumerate(datelist)])
      self.number_of_time_records   = len(lst)
      if self.number_of_time_records==0:
         return

      self.datetimes = []
      self.filelist  = []
      self.steplist  = []
      self.objects  = []
      for dto,i in lst:
         f     = filelist[i]
         step  = steplist[i]
         self.datetimes.append(dto)
         self.filelist.append(f)
         self.steplist.append(step)

         # get binary file object for each file
         print( 'Getting object for step %d (%d - %d)...' %(step,np.min(steplist),np.max(steplist)) )
         nbi   = nextsim_binary_info(f,**kwargs)
         self.objects.append(nbi)

         if len(self.objects)==1:
            # gmsh boundary object applies to all (boundary nodes don't move),
            # and takes some time to create,
            # so just make it once
            self.boundary  = nbi.mesh_info.boundary
            if self.boundary is not None:
               kwargs.update({'boundary':self.boundary})

            # reading log file takes a bit of time too
            self.options   = nbi.options
            if self.options is not None:
               kwargs.update({'options':self.options})
      # ======================================================================

      self.variables = self.objects[0].variables
      t1 = time.time()
      print('\nSorting/initialisation time: %f seconds.'%(t1-t0))
      return
   # =========================================================================


   # =========================================================================
   def nearestDate(self, pivot):
      """
      dto,time_index = self.nearestDate(dto0)
      dto0  = datetime.datetime objects
      dto   = datetime.datetime objects - nearest value in self.datetimes to dto0
      time_index: dto=self.datetimes[time_index]
      """
      dto         = min(self.datetimes, key=lambda x: abs(x - pivot))
      time_index  = self.datetimes.index(dto)
      return dto,time_index
   # =========================================================================


   # =========================================================================
   def get_var(self,varname,time_index=0):
      return self.objects[time_index].get_var(varname,**kwargs)
   # =========================================================================


   # =============================================================================
   def plot_mesh(self,time_index=0,**kwargs):
      return self.objects[time_index].plot_mesh(**kwargs)
   # =============================================================================


   # =============================================================================
   def plot_var(self,vname,time_index=0,**kwargs):
      nbi   = self.objects[time_index]
      return nbi.plot_var(vname,**kwargs)
   # =============================================================================


   # =============================================================================
   def plot_var_all(self,vname,figdir='.',**kwargs):

      import os
      if 'pobj' in kwargs:
         pobj  = kwargs['pobj']
      else:
         pobj  = None
      
      if not os.path.exists(figdir):
         os.mkdir(figdir)
      figdir2  = figdir+'/'+vname
      if not os.path.exists(figdir2):
         os.mkdir(figdir2)

      for i,nbi in enumerate(self.objects):

         kwargs.update({'pobj':pobj})
         kwargs.update({'show':False})
         figname  = figdir2+'/'+nbi.basename+'.png'

         pobj  = nbi.plot_var(vname,figname=figname,**kwargs)
         pobj.ax.cla()

      from matplotlib import pyplot as plt
      plt.close(pobj.fig)
      return
   # =============================================================================
