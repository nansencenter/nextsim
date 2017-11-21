# ===============================================================
def xyz_to_lonlat(x,y,z):
   import numpy as np
   PI       = np.pi
   radius   = np.sqrt(pow(x,2.)+pow(y,2.)+pow(z,2.));
   lat      = np.arcsin(z/radius)*(180./PI);
   lon      = np.arctan2(y,x);
   lon      = lon-2*PI*np.floor(lon/(2*PI));
   lon      = lon*(180./PI);
   return lon,lat
# ===============================================================


# ===============================================================
def jacobian(x0,y0,x1,y1,x2,y2):
    return (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0)
# ===============================================================


# ===============================================================
def get_areas(mesh_obj):

   import numpy as np
   nodes_x,nodes_y   = mesh_obj.get_nodes_xy()
   indices           = mesh_obj.get_indices(eltype="triangles",numbering='gmsh',asvector=False)
   areas             = []

   for inds in indices:
      x  = []
      y  = []
      for n in inds:
         x.append(nodes_x[n])
         y.append(nodes_y[n])

      areas.append( .5*np.abs(jacobian(x[0],y[0],x[1],y[1],x[2],y[2])) )

   return areas
# ===============================================================


# ===============================================================
def get_resolution(mesh_obj):

   import numpy as np
   areas = get_areas(mesh_obj)
   return np.sqrt(np.mean(areas))
# ===============================================================


# ===============================================================
def default_pyproj():
   import pyproj
   import numpy as np

   ecc      = 0.081816153
   a        = 6378.273e3 #equatorial radius in m
   b        = a*np.sqrt(1-pow(ecc,2))
   lat0     = 90.
   lon0     = -45.
   lat_ts   = 60.
   return pyproj.Proj(proj='stere',a=a,b=b,\
                  lon_0=lon0,lat_0=lat0,lat_ts=lat_ts)
# ===============================================================


# ===============================================================
def mppfile_to_pyproj(mppfile=None):

   if mppfile is None:
      return default_pyproj()

   mf       = open(mppfile)
   lines    = mf.readlines()
   mf.close()

   import pyproj
   import numpy as np

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
# ===============================================================


# ===============================================================
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
# ===============================================================


# ===============================================================
def get_array(vname,binfile,vnames,vtypes):

   if vname not in vnames:
      raise ValueError(vname+ ' not in '+binfile)

   import struct
   import numpy as np

   sizes = {'f':4,'i':4,'d':8}

   f  = open(binfile,'rb')
   for v in vnames:
      # print(v)
      fmt   = vtypes[v]
      size  = sizes[fmt]

      # determine record length
      N  = f.read(4)
      N  = struct.unpack('i',N)[0]

      if v!=vname:
         f.seek(N*size,1)
      else:
         data  = f.read(N*size)
         data  = np.array(struct.unpack(N*fmt,data))
         break

   f.close()
   return data
# ===============================================================


# ===============================================================
def get_variable_lengths(binfile,vnames,vtypes):

   import struct
   sizes = {'f':4,'i':4,'d':8}

   vlens = {}
   f  = open(binfile,'rb')
   for v in vnames:
      # print(i,v)
      fmt   = vtypes[v]
      size  = sizes[fmt]

      N  = f.read(4)
      N  = struct.unpack('i',N)[0]
      vlens.update({v:N})
      f.seek(N*size,1)

   f.close()
   return vlens
# ===============================================================


# ===============================================================
def get_arrays(binfile,vnames,vtypes):

   import struct
   import numpy as np

   sizes = {'f':4,'i':4,'d':8}
   out   = {}

   f  = open(binfile,'rb')
   for v in vnames:
      # print(v)
      fmt   = vtypes[v]
      size  = sizes[fmt]

      # determine record length
      N  = f.read(4)
      N  = struct.unpack('i',N)[0]

      data  = f.read(N*size)
      out.update({v:np.array(struct.unpack(N*fmt,data))})

   f.close()
   return out
# ===============================================================


# ===============================================================
def define_grid(mesh_obj,resolution=None):

   # ============================================================
   # initialise
   if resolution is None:
      if hasattr(mesh_obj,'resolution'):
         res   = .5*mesh_obj.resolution
      else:
         raise ValueError('please specify resolution as mesh_obj does not have a resolution attribute')
   else:
      res   = float(resolution)

   import numpy as np

   xmin  = mesh_obj.xmin
   xmax  = mesh_obj.xmax
   ymin  = mesh_obj.ymin
   ymax  = mesh_obj.ymax
   xav   = .5*(xmin+xmax)
   yav   = .5*(ymin+ymax)
   
   # increase Dx,Dy to give correct resolution
   Dx    = xmax-xmin
   Dy    = ymax-ymin
   nx    = int(np.ceil(Dx/res))
   ny    = int(np.ceil(Dy/res))
   Dx    = nx*res
   Dy    = ny*res
   xmin  = xav-.5*Dx
   xmax  = xav+.5*Dx
   ymin  = yav-.5*Dy
   ymax  = yav+.5*Dy
   # ============================================================

   grid_params = {}
   grid_params.update({'xmax':xmax})
   grid_params.update({'ymax':ymax})
   grid_params.update({'xmin':xmin})
   grid_params.update({'ymin':ymin})
   grid_params.update({'nx':nx})
   grid_params.update({'ny':ny})

   return grid_params
# ===============================================================


# ===============================================================
def interp_mesh_to_points(mesh_obj,xout,yout,data_in,**kwargs):
      
   import bamg # wrapper for some functions in bamg c++ library
   import numpy as np

   nodes_x,nodes_y   = mesh_obj.get_nodes_xy(dtype='double')
   indices           = mesh_obj.get_indices(eltype="triangles",numbering='bamg',asvector=True)

   # make sure input data is list of arrays of type double
   data     = []
   dtype    = type(data_in)
   islist   = dtype==type([])
   isdict   = dtype==type({})

   shp   = xout.shape
   sz    = xout.size

   # make xout,yout into vectors
   X     = np.reshape(xout.astype('double'),(sz))
   Y     = np.reshape(yout.astype('double'),(sz))

   if isdict:
      for key in data_in:
         data.append( np.array(data_in[key],dtype='double') )
   else:
      for dat in data_in:
         data.append( np.array(dat,dtype='double') )

   print('calling bamg interp')
   data_out = bamg.interpMeshToPoints(indices,nodes_x,nodes_y,
                                       data,X,Y,**kwargs)

   if isdict:
      # return as dict of arrays the same shape as xout,yout
      out   = {}
      for i,key in enumerate(data_in):
         out.update({ key:np.reshape(data_out[i],shp) })
      return out
   else: 
      out   = []
      for dat in data_out:
         # return as list of arrays the same shape as xout,yout
         out.append( np.reshape(dat,shp) )
      return out
# ===============================================================


# ===============================================================
def sort_program_options(po0):
   po = {}
   for key0 in po0:
      val         = po0[key0]
      key,subkey  = key0.split('.')
      if key not in po:
         po.update({key:{}})
      po[key].update({subkey:val})
   return po
# ===============================================================


# ===============================================================
def read_nextsim_log(logfile):

   if logfile is None:
      return

   from datetime import datetime as DT
   fid   = open(logfile,'r')
   lines = fid.readlines()
   fid.close()

   options        = {}
   keys           = []
   config_files   = []
   badlines       = []
   for lin in lines:
      # print(lin)

      if '#---' in lin:
         # header of section
         # - these will be the keys of a dictionary
         key   = lin.split('-')[-1]
         key   = key.replace(' ','_')
         key   = key.split()[0].lower()# removes \n and makes lower case
         keys.append(key) 
         options.update({key:{}})
         continue

      elif ']=' in lin:
         # multiple config files
         config_files.append(lin.split(']=')[1])
         options.update({'config_files':config_files})
         continue

      else:
         ls = lin.split()

         # ======================================================
         if len(ls)==2:
            # normal situation
            key,val0 = ls

         elif len(ls)==1:
            # key name is probably too long - skip
            badlines.append(lin)
            continue

         elif ls[0].lower()=='build':
            # special case
            key      = 'build_date'
            keywidth = lin.find(ls[2])
            fmt      = '%Y-%b-%d %H:%M:%S'
            ds       = ls[2]+' '+ls[3]
            val0     = DT.strptime(ds,fmt)

         else:
            key      = lin[:keywidth].replace(' ','_').strip('_')
            val0     = lin[keywidth:].replace(' ','_')
         # ======================================================


         # ======================================================
         try:
            # convert to numeric value
            val   = float(val0)
         except:
            # keep as string
            val   = val0
         # ======================================================

         options[keys[-1]].update({key:val})

   options['program_options'] = sort_program_options(options['program_options'])

   if len(badlines)>0:
      print("WARNING: couldn't parse these lines (either too long or no value):\n")
      for lin in badlines:
         print(lin)
      print("\n")

   return options
# ===============================================================
