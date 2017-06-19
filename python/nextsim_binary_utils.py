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
      self.nodes_lon    = []
      self.nodes_lat    = []
      for n in range(self.num_nodes):
         lin            = fid.readline()
         ident,c1,c2,c3 = lin.split()
         self.nodes_id.append(int(ident))

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
            verts.append(int(slin[jcpt]))
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
