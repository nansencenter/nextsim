## make_nextwim_grid.py
## Author: Timothy Williams
## Date: 20160916
import numpy as np

# ===============================================================
def sav2bin(aid,arr,fields,nbytes=8):
   import struct

   key   = arr.keys()[0]
   fields.append(key)

   X     = arr[key]
   nX    = X.size

   print("Range of "+key,X.min(),X.max())

   if nbytes==8:
      ss = nX*'d'
   else:
      ss = nX*'f'

   # reshape to fortran order
   A  = X.reshape((nX),order='F')
   aid.write(struct.pack(ss,*A))
   
   return
# ===============================================================

# ===============================================================
def write_bfile(bid,fields,nx,ny,nbytes=8):

   nrecs = len(fields)

   bid.write("%2.2i       Nrecs    # Number of records\n" %(nrecs))
   bid.write("01       Norder   # Storage order [column-major (F/matlab) = 1, row-major (C) = 0]\n")
   bid.write("%4.4i     nx       # Record length in x direction (elements)\n" %(nx))
   bid.write("%4.4i     ny       # Record length in y direction (elements)\n" %(ny))
   bid.write("%2.2i       nbytes       # bytes per entry\n" %(nbytes))
   bid.write("\n")
   bid.write("Record number and name:\n")

   """
   Rest of file:
   01       qlon
   02       qlat
   03       plon
   04       plat
   05       ulon
   06       ulat
   07       vlon
   08       vlat
   09       scuy
   10       scvx
   11       LANDMASK
   """

   for recno,vname in enumerate(fields):
      ss = "%2.2i       %s\n" %(recno+1,vname)
      bid.write(ss)

   return
# ===============================================================


# =======================================================================
def get_depth_mask(plon,plat,ncfil=None,**kwargs):

   if ncfil is None:
      ncfil = 'ETOPO/ETOPO_Arctic_10arcmin.nc' # ~20km
      # ncfil = 'ETOPO/ETOPO_Arctic_5arcmin.nc'  # ~10km
      # ncfil = 'ETOPO/ETOPO_Arctic_2arcmin.nc'  # ~4km
      # ncfil = 'ETOPO/ETOPO_Arctic_1arcmin.nc'  # ~2km

   #interp depth to plon,plat (centres of wim grid cells)
   import mod_reading as mr
   nci      = mr.nc_getinfo(ncfil)
   depth    = nci.interp2points('z', (plon,plat),**kwargs) # masked array
   landmask = np.ones(depth.shape)
   good     = np.logical_not(depth.mask)
   data     = depth.data[good]
   lmgood   = landmask[good] 

   lmgood[data<0.]   = 0. # z<0 is water
   landmask[good]    = lmgood

   return depth,landmask
# =======================================================================


# =======================================================================
def get_connected(landmask,TEST_CONN=True):
   from skimage import measure as msr
   from skimage import morphology as morph
   
   # find connected regions
   # - land (landmask==1) is the background
   # - connectivity=1, so diagonal connectivity doesn't count
   labels,Nlabels = msr.label(landmask,background=1,return_num=True,connectivity=1)

   # ==================================================
   # find largest connected region
   maxlen   = 0
   maxlab   = 0
   for n in range(1,Nlabels):
      N  = len(labels[labels==n])
      # print(n,N)
      if N>maxlen:
         maxlen   = N
         maxlab   = n
         # print(n,N)

   Landmask = 1+0*landmask
   Landmask[labels==maxlab]   = 0.
   # ==================================================


   # ==================================================
   # open to avoid "singular points"
   selem = np.ones((3,3))
   if 1:
      selem=None # default (corners=0)
   Landmask = morph.binary_opening(Landmask.astype('int'),selem=selem).astype(int)

   if 0:
      # plot original landmask
      fig1  = plt.figure()
      ax1   = fig1.add_subplot(1,1,1)
      I1    = ax1.imshow(landmask.transpose(),origin='lower')
      fig1.colorbar(I1)
      fig1.show()

   if TEST_CONN:
      # plot processed landmask
      fig3  = plt.figure()
      ax3   = fig3.add_subplot(1,1,1)
      I3    = ax3.imshow(Landmask.transpose(),origin='lower')
      fig3.colorbar(I3)
      plt.show(fig3)

   if 0:
      # plot all labels found
      fig2  = plt.figure()
      ax2   = fig2.add_subplot(1,1,1)
      I2 = ax2.imshow(labels,cmap='spectral')
      fig2.colorbar(I2)
      plt.show(fig2)

   return Landmask
# =======================================================================

class wim_grid:
   def __init__(self,\
                  gmsh_mesh=None,\
                  gmsh_mesh_file=None,\
                  grid_info_file=None,\
                  grid_params=None,\
                  mapping=None,\
                  mppfile=None):
      
      self.gmsh_mesh       = None
      self.gmsh_mesh_file  = None
      self.grid_info_file  = None

      self.use_gmsh        = False
      selg.mapping         = mapping

      # ====================================================================================
      # get gmsh_mesh if wanted
      if gmsh_mesh is not None:
         self.use_gmsh  = True
         self.gmsh_mesh = gmsh_mesh
      elif gmsh_mesh_file is not None:
         import nextsim_binary_utils as nbu
         self.use_gmsh    = True
         self.gmsh_mesh   = nbu.gmsh_mesh(gmsh_mesh_file,mapping=mapping,mppfile=mppfile)
      # ====================================================================================


      # ====================================================================================
      # get grid_params
      if self.use_gmsh:
         # stereographic projection already defined
         self.mapping   = gmsh_mesh.mapping
         grid_params    = gmsh_mesh.export_wim_grid_params()
      else:
         # may still need to set projection to use
         if self.mapping is None:
            self.mapping   = nbu.mppfile_to_pyproj(mppfile=mppfile)

         # check for grid_params
         if grid_info_file is not None:
            gridname,grid_params = self.read_grid_info_file(grid_info_file)
         elif grid_params is None:
            raise ValueError('Please pass in one of gmsh_mesh, gmsh_mesh_file, grid_info_file or grid_params')
      # ====================================================================================


      # ====================================================================================
      # should now have grid_params
      # - make the grid
      self.make_grid(grid_params)
      # ====================================================================================

      return
      
      # # - set the land mask
      # self.set_land_mask()

      # # - write to binary file
      # self.write_binary(gridname)
      # # ====================================================================================

      # return
   # =======================================================================

   # =======================================================================
   def read_grid_info_file(filename):

      # read the file
      f     = open(filename)
      lines = f.readlines()
      f.close()
      res         = float(lines[1].split()[0])
      gridname    = lines[0].split()[0]
      grid_params = {}

      vertices = []
      for lin in lines[2:]:
         lon   = float(lin.split()[0])
         lat   = float(lin.split()[1])
         vertices.append((lon,lat))
      
      for lon,lat in vertices:
         x,y   = self.mapping(lon,lat)
         xc.append(x)
         yc.append(y)

      xmin  = np.min(xc)
      xmax  = np.max(xc)
      ymin  = np.min(yc)
      ymax  = np.max(yc)
      Dx    = xmax-xmin
      Dy    = ymax-ymin

      # reduce xmax,ymax to give correct resolution
      nx    = int(np.floor(Dx/res))
      ny    = int(np.floor(Dy/res))
      if 0:
         xav   = (xmin+xmax)/2.
         yav   = (ymin+ymax)/2.
         xmin  = xav-nx*res/2.
         ymin  = yav-ny*res/2.

      # final output
      grid_params.update({'xmax':xmin+nx*res})
      grid_params.update({'ymax':ymin+ny*res})
      grid_params.update({'xmin':xmin})
      grid_params.update({'xmin':xmin})
      grid_params.update({'nx':ny})
      grid_params.update({'ny':ny})

      print('Resolution (km): '+str(res/1.e3))
      print('x range (km)   : '+str((xmax-xmin)/1.e3))
      print('y range (km)   : '+str((ymax-ymin)/1.e3))
      print('nx             : '+str(nx))
      print('ny             : '+str(ny))

      return gridname,grid_params
   # =======================================================================

   # =======================================================================
   def make_grid(self):

      nx    = self.grid_params['nx']
      ny    = self.grid_params['ny']
      xmin  = self.grid_params['xmin']
      ymin  = self.grid_params['ymin']
      xmax  = self.grid_params['xmax']
      ymax  = self.grid_params['ymax']

      # corners: (nx+1)*(ny+1)
      # X increases down rows, Y increases along columns
      qx    = np.linspace(xmin,xmax,nx+1)
      qy    = np.linspace(ymin,ymax,ny+1)
      qY,qX = np.meshgrid(qy,qx)
      self.arrays.update({'qx':qX})
      self.arrays.update({'qy':qY})

      # centres: nx*ny
      px    = .5*(qx[1:]+qx[:-1])
      py    = .5*(qy[1:]+qy[:-1])
      pY,pX = np.meshgrid(py,px)
      self.arrays.update({'px':pX})
      self.arrays.update({'py':pY})

      # side points: (nx+1)*ny
      uY,uX = np.meshgrid(py,qx)
      self.arrays.update({'ux':uX})
      self.arrays.update({'uy':uY})

      # up/down points: nx*(ny+1)
      vY,vX = np.meshgrid(qy,px)
      self.arrays.update({'vx':vX})
      self.arrays.update({'vy':vY})

      # ================================================
      # grid cell sizes
      self.arrays.update({'scuy':np.zeros((nx+1,ny))+res})
      self.arrays.update({'scvx':np.zeros((nx,ny+1))+res})
      self.arrays.update({'scpy':np.zeros((nx,ny))+res})
      self.arrays.update({'scpx':np.zeros((nx,ny))+res})

      self.arrays.update({'scuy':scuy})
      self.arrays.update({'scvx':scvx})
      self.arrays.update({'scp2':scp2})
      # ================================================


      return
   # =======================================================================
   # =======================================================================
   def write_binary(gridname,outdir='.'):

      # ===============================================================
      # write to binary
      aid      = open(outdir+'/'+gridname+'.a','wb')
      fields   = [] # list of variable names

      # corners: (nx+1)*(ny+1)
      # X increases down rows, Y increases along columns
      if not self.USE_LL:
         sav2bin(aid,{'qx':self.arrays['qx']},fields)
         sav2bin(aid,{'qy':self.arrays['qy']},fields)
      else:
         qlon,qlat   = self.mapping(self.arrays['qx'],self.arrays['qy'],inverse=True)
         sav2bin(aid,{'qlon':qlon},fields)
         sav2bin(aid,{'qlat':qlat},fields)

      # centres: nx*ny
      px    = .5*(qx[1:]+qx[:-1])
      py    = .5*(qy[1:]+qy[:-1])
      pY,pX = np.meshgrid(py,px)
      if not self.USE_LL:
         sav2bin(aid,{'px':pX},fields)
         sav2bin(aid,{'py':pY},fields)
      else:
         plon,plat   = self.mapping(pX,pY,inverse=True)
         sav2bin(aid,{'plon':plon},fields)
         sav2bin(aid,{'plat':plat},fields)

      # side points: (nx+1)*ny
      if not self.USE_LL:
         sav2bin(aid,{'ux':self.arrays['ux']},fields)
         sav2bin(aid,{'uy':self.arrays['uy']},fields)
      else:
         ulon,ulat   = self.mapping(self.arrays['ux'],self.arrays['uy'],inverse=True)
         sav2bin(aid,{'ulon':ulon},fields)
         sav2bin(aid,{'ulat':ulat},fields)

      # up/down points: nx*(ny+1)
      if not USE_LL:
         sav2bin(aid,{'vx':self.arrays['vx']},fields)
         sav2bin(aid,{'vy':self.arrays['vy']},fields)
      else:
         vlon,vlat   = self.mapping(self.arrays['vx'],self.arrays['vy'],inverse=True)
         sav2bin(aid,{'vlon':vlon},fields)
         sav2bin(aid,{'vlat':vlat},fields)

      # ================================================
      # grid cell sizes

      if self.USE_LL:
         # CHANGE TO DISTANCES ON THE ELLIPSOID
         scuy  = np.zeros((nx+1,ny))
         scvx  = np.zeros((nx,ny+1))
         scpy  = np.zeros((nx,ny))
         scpx  = np.zeros((nx,ny))

         # uy,py:
         # y increases in j dirn, so get distance between cols
         for j in range(ny):
            # loop over cols
            scuy[:,j]   = geod.inv(qlon[:,j+1],qlat[:,j+1],qlon[:,j],qlat[:,j])[2]
            scpy[:,j]   = geod.inv(vlon[:,j+1],vlat[:,j+1],vlon[:,j],vlat[:,j])[2]

         # vx:
         # x increases in i dirn, so get distance between rows
         for i in range(nx):
            # loop over rows
            scvx[i,:]   = geod.inv(qlon[i+1,:],qlat[i+1,:],qlon[i,:],qlat[i,:])[2]
            scpx[i,:]   = geod.inv(ulon[i+1,:],ulat[i+1,:],ulon[i,:],ulat[i,:])[2]

         scp2  = scpx*scpy
         sav2bin(aid,{'scuy':scuy},fields)
         sav2bin(aid,{'scvx':scvx},fields)
         sav2bin(aid,{'scp2':scp2},fields)
      # ================================================

      sav2bin(aid,{'LANDMASK':landmask},fields)
      aid.close()

      bid   = open(outdir+'/'+gridname+'.b','w')
      write_bfile(bid,fields,nx,ny)
      bid.close()
      print("\n grid-files "+outdir+'/'+gridname+'.[a,b] saved\n')

      return
   # =======================================================================


   # =======================================================================
   def set_landmask(self,depth_file=None,select_opt=1):
      if self.use_gmsh:
         if select_opt==1:
            wet   = self.gmsh_mesh.boundary.iswet(self.arrays['px'],self.arrays['py'])
         else:
            qwet  = self.gmsh_mesh.boundary.iswet(self.arrays['qx'],self.arrays['qy'])
            wet   = np.logical_and(qmask[:nx-1,1:],qmask[:nx-1,:ny-1])
            wet   = np.logical_and(wet,qmask[1:,1:])
            wet   = np.logical_and(wet,qmask[1:,:ny-1])

         self.arrays['landmask'] = np.logical_not(wet)
      elif depth_file is not None:
         plon,plat   = self.mapping(self.arrays['px'],self.arrays['py'])
         land_mask   = get_depth_mask(plon,plat,ncfil=depth_file,mapping=self.mapping)

         # make sure it's connected
         self.arrays['landmask'] = get_connected(land_mask,TEST_CONN=False)

      return
   # =======================================================================
