import os,sys
import numpy as np

def xyz_to_lonlat(x,y,z):
   PI       = np.pi
   radius   = np.sqrt(pow(x,2.)+pow(y,2.)+pow(z,2.));
   lat      = np.arcsin(z/radius)*(180./PI);
   lon      = np.arctan2(y,x);
   lon      = lon-2*PI*np.floor(lon/(2*PI));
   lon      = lon*(180./PI);
   return lon,lat

def default_pyproj():
   import pyproj
   ecc      = 0.081816153
   a        = 6378.273e3 #equatorial radius in m
   b        = a*np.sqrt(1-pow(ecc,2))
   lat0     = 90.
   lon0     = -45.
   lat_ts   = 60.
   return pyproj.Proj(proj='stere',a=a,b=b,\
                  lon_0=lon0,lat_0=lat0,lat_ts=lat_ts)

def mppfile_to_pyproj(mppfile=None):

   if mppfile is None:
      return default_pyproj()

   mf       = open(mppfile)
   lines    = mf.readlines()
   mf.close()

   import pyproj

   # shape of earth
   ecc   = float(lines[-1].split()[0])
   a     = 1e3*float(lines[-2].split()[0]) # convert to m
   b     = a*np.sqrt(1-pow(ecc,2))

   # stere info
   lat0,lon0,lat_ts  = lines[1].split()[:3]
   lon0              = float(lon0)
   lat0              = float(lat0)
   lat_ts            = float(lat_ts)

   rotation = float(lines[2].split()[0])

   print(a,b,ecc)
   print(rotation,lat0,lat_ts)
   
   return pyproj.Proj(proj='stere',a=a,b=b,\
                  lon_0=rotation,lat_0=lat0,lat_ts=lat_ts)


def read_file_info(file_info):
      f  = open(file_info,'r')
      lines = f.readlines()
      f.close()

      variables       = []
      variable_types  = {} # list of variable types: 'f' or 'd' (single or double) or 'i' (integer)
      record_numbers  = {}

      for recno,lin in enumerate(lines):
         ss = lin.split()
         variables.append(ss[0])
         record_numbers.update({ss[0]:recno})
         if ss[1]=='double' or ss[1]=='8':
            variable_types.update({ss[0]:'d'})
         elif ss[1]=='single' or ss[1]=='float' or ss[1]=='4':
            variable_types.update({ss[0]:'f'})
         elif ss[1]=='int':
            variable_types.update({ss[0]:'i'}) # signed integer
         else:
            raise ValueError('unknown data type in '+file_info+': '+ss[1])

      return variables,record_numbers,variable_types

class nextsim_mesh_info:
   def __init__(self,mesh_file):

      self.mesh_file       = mesh_file
      self.mesh_file       = os.path.abspath(mesh_file)
      self.mesh_file_info  = self.mesh_file.replace('.bin','.dat')

      self.variables,self.record_numbers,self.variable_types  = read_file_info(self.mesh_file_info)

      return

class nextsim_binary_info:
   def __init__(self,data_file):

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
      self.mesh_info = nextsim_mesh_info(mesh_file)

      #TODO read nextsim.log (mppfile to initialise projection, original mesh_file)
      #TODO or pass them in as arguments?

      return

   def _read_field_info(self):

      import struct
      from datetime import datetime as DT,timedelta as TD

      # get the variables and their record no
      self.variables,self.record_numbers,self.variable_types  = read_file_info(self.data_file_info)

      # Time format
      fmt   = self.variable_types['Time']
      if fmt == 'f':
         size = 4
      elif fmt == 'd':
         size = 8

      f  = open(self.data_file,'rb')
      f.read(4)
      Time  = f.read(size)
      f.close()

      Time  = struct.unpack(fmt,Time)[0] #days
      
      refdate        = DT(1900,1,1)
      self.datetime  = refdate+TD(Time)
      self.datetimes = [self.datetime]

      return

class mesh_physical_name:
   def __init__(self,name,ident,topodim):
      self.name      = name
      self.ident     = ident
      self.topodim   = topodim
      return

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

class gmsh_boundary:
   def __init__(self,exterior,islands=None,open_boundaries=None,coastal_boundaries=None):
      self.exterior_polygon   = exterior           # Shapely Polygon
      self.island_polygons    = islands            # List of shapely Polygons
      self.open_boundaries    = open_boundaries    # List of shapely LineStrings
      self.coastal_boundaries = coastal_boundaries # List of shapely LineStrings
      return

   def plot(self,sort=True,show=True):
      from matplotlib import pyplot as plt
      fig   = plt.figure()
      ax    = fig.add_subplot(1,1,1)
      sfac  = 1e-3 #m to km
      lw    = 1.35

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

      return fig,ax

class gmsh_mesh:

   # ================================================================
   def __init__(self,meshfile,ordering='GMSH',mppfile=None,mapping=None):
      # TODO: add BAMG ordering permutation

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
         self.mapping   = mppfile_to_pyproj(self.mppfile)

      self.stereographic_projection()

      self.get_boundary()
      return
   # ================================================================

   # ================================================================
   def read_meshfile(self):

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

         lon,lat  = xyz_to_lonlat(float(c1),float(c2),float(c3))
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


   # ================================================================
   def get_indices(self,eltype="triangles"):

      indices  = []
      if eltype=="triangles":
         els   = self.triangles
      else:
         els   = self.edges
      for el in els:
         indices.append(el.vertices)
         
      return np.array(indices)
   # ================================================================


   # ================================================================
   def plot_mesh(self,plotter="pyplot",**kwargs):

      if plotter=="pyplot":
         # patches used in self.plot_data() is quicker than trying to
         # plot individual triangles
         self.plot_data(plot_grid=True,plot_coast=True,**kwargs)

      elif plotter=="mlab":
         from mayavi import mlab
         z        = 0*self.nodes_x
         indices  = self.get_indices("triangles")
         mlab.triangular_mesh(self.nodes_x,self.nodes_y,z,indices,\
               representation='wireframe',color=(0,0,0))
         mlab.view(0,0)
         mlab.show()

      return
   # ================================================================


   # ================================================================
   def plot_data(self,data=None,clabel=None,plot_grid=False,plot_coast=True,show=True,colorbar=True):
      from matplotlib import patches,cm,collections
      from matplotlib import pyplot as plt
      fig   = plt.figure()
      ax    = fig.add_subplot(1,1,1)
      sfac  = 1e-3 #m to km
      ax.set_xlim([sfac*self.nodes_x.min(),sfac*self.nodes_x.max()])
      ax.set_ylim([sfac*self.nodes_y.min(),sfac*self.nodes_y.max()])

      indices  = self.get_indices("triangles")
      if data is None:
         # test data: plot latitude
         data     = np.mean(self.nodes_lat[indices],axis=1)
         clabel   = 'Latitude'

      patch_list  = []
      for triangle in self.triangles:
         ccl   = triangle.get_coords(sfac*self.nodes_x,sfac*self.nodes_y)
         ccl.append(ccl[0]) # close the contour
         patch_list.append(patches.Polygon(ccl,True))

      pc = collections.PatchCollection(patch_list, cmap=cm.jet, alpha=.5)
      pc.set_array(data)

      if plot_grid:
         pc.set_edgecolor('k')

      ax.add_collection(pc)
      ax.set_xlabel('x, km')
      ax.set_ylabel('y, km')

      if colorbar:
         cbar  = fig.colorbar(pc)
         if clabel is not None:
            cbar.set_label(clabel,rotation=270)

      if plot_coast:
         self.plot_boundary(fig=fig,ax=ax,sort=True,show=False)

      if show:
         fig.show()
      
      return fig,ax
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
      lens  = []
      polys = []
      for ll in lm:
         pp = shg.Polygon(ll)
         polys.append(pp)
         lens.append(len(pp.exterior.coords))

      if len(polys)==1:
         islands  = None
      else:
         lst      = sorted([(e,i) for i,e in enumerate(lens)],reverse=True)
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
   def plot_boundary(self,sort=False,show=True,fig=None,ax=None):
      from matplotlib import collections
      from matplotlib import pyplot as plt

      sfac     = 1e-3 # m to km
      haveax   = True
      if fig is None:
         fig      = plt.figure()
         ax       = fig.add_subplot(1,1,1)
         haveax   = False
      elif ax is None:
         ax       = fig.add_subplot(1,1,1)
         haveax   = False

      if not haveax:
         # if ax was passed in, assume these are set outside
         ax.set_xlim([sfac*self.nodes_x.min(),sfac*self.nodes_x.max()])
         ax.set_ylim([sfac*self.nodes_y.min(),sfac*self.nodes_y.max()])
         ax.set_xlabel('x, km')
         ax.set_ylabel('y, km')

      if sort:
         if len(self.PhysicalNames)==0:
            plabs = []
            cpt   = 0
            while (self.edges[cpt].physical not in plabs) and (len(plabs)<2):
               plabs.append(self.edges[cpt].physical)
               cpt  += 1
            coast_val   = plabs[0]
            open_val    = plabs[1]
         else:
            for ident in self.PhysicalNames:
               if self.PhysicalNames[ident].name=="open":
                  open_val = ident
               if self.PhysicalNames[ident].name=="coast":
                  coast_val   = ident

      if not sort:
         line_list  = []
         for edge in self.edges:
            ccl   = edge.get_coords(sfac*self.nodes_x,sfac*self.nodes_y)
            line_list.append(ccl)
         lcols = [collections.LineCollection(line_list,colors='k')]
      else:
         open_list   = []
         coast_list  = []
         for edge in self.edges:
            ccl   = edge.get_coords(sfac*self.nodes_x,sfac*self.nodes_y)
            if edge.physical==open_val:
               open_list.append(ccl)
            else:
               coast_list.append(ccl)

         lcols = [collections.LineCollection(open_list ,colors='c'),\
                  collections.LineCollection(coast_list,colors='k')]

      for lc in lcols:
         lc.set_linewidth(1.)
         ax.add_collection(lc)

      if show:
         fig.show()
      
      return fig,ax
   # ================================================================
