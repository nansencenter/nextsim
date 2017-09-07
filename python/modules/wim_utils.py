import numpy as np
import os,sys
import struct
import fns_plot_data as Fplt
import datetime as dtm

class plot_object:
   def __init__(self):
      from matplotlib import pyplot as plt
      self.fig   = plt.figure()
      self.ax    = self.fig.add_subplot(1,1,1)
      return

##############################################################
def fn_bfile_info(bfile):
   # routine to get output fields from binary files:

   ###########################################################
   # get info like dimensions and variable names from .b file
   bid   = open(bfile,'r')
   lines = bid.readlines()
   bid.close()

   int_list = ['nx','ny','Nrecs','Norder'] # these are integers not floats
   binfo    = {'recnos':{}}  # dictionary with info about fields in corresponding .a file

   do_vlist = 0
   count    = 0
   for lin in lines:
      ls = lin.split()
      if ls != []:
         # skip blank lines

         # if not at "Record number and name"
         if not(ls[0]=='Record' and ls[1]=='number'):
            val   = ls[0]
            key   = ls[1]
            if not do_vlist:
               if key in int_list:
                  val   = int(val)
               elif ("T" in val) and ("Z" in val):
                  import datetime as dtm
                  val   = dtm.datetime.strptime(val,"%Y%m%dT%H%M%SZ")
               else:
                  val   = float(val)
               binfo.update({key:val})
            else:
               binfo['recnos'].update({key:count})
               count   += 1
         else:
            # have got down to list of variables
            do_vlist = 1

   if binfo['Nrecs']!=len(binfo['recnos']):
      raise ValueError('Inconsistent number of records in file: '+bfile)

   return binfo
##############################################################


##############################################################
def get_array(fid,recno,nx,ny,fmt_size=4,order='F'):
   # routine to get the array from the .a (binary) file
   # * fmt_size = size in bytes of each entry)
   #   > default = 4 (real*4/single precision)
   recs     = nx*ny
   rec_size = recs*fmt_size
   fid.seek(recno*rec_size,0) # relative to start of file
   #
   if fmt_size==4:
      fmt_py   = 'f' # python string for single
   else:
      fmt_py   = 'd' # python string for double

   data  = fid.read(rec_size) # no of bytes to read
   fld   = struct.unpack(recs*fmt_py,data)
   fld   = np.array(fld)
   fld   = fld.reshape((nx,ny),order=order)

   return fld
##############################################################


###########################################################
def check_names(vname,variables,stop=True):

   if vname in variables:
      return vname

   lists = []

   # ice conc alt names
   List   = ['ficem','fice','ice_conc','icec','cice','area',\
                  'concentration','sea_ice_concentration']
   if vname in List:
      for v in variables:
         if v in List:
            return v

   # ice thick alt names
   List  = ['hicem','hice','ice_thick','icetk','iceh',\
            'sea_ice_thickness','thickness']
   if vname in List:
      for v in variables:
         if v in List:
            return v

   # floe size alt names
   List  = ['dfloe','dmax','Dfloe','Dmax']
   if vname in List:
      for v in variables:
         if v in List:
            return v

   # wave stress: x component
   List  = ['taux','tau_x','taux_waves']
   if vname in List:
      for v in variables:
         if v in List:
            return v

   # wave stress: y component
   List  = ['tauy','tau_y','tauy_waves']
   if vname in List:
      for v in variables:
         if v in List:
            return v

   # swh
   List  = ['Hs','hs','swh','significant_wave_height']
   if vname in List:
      for v in variables:
         if v in List:
            return v

   if stop:
      print('Failed to get '+vname)
      print('\nAvailable variables:')
      for v in variables:
         print(v)

      print('')
      raise ValueError(vname+' not in variable list')
   else:
      return ''
###########################################################


##############################################################
def fn_read_general_binary(afile,vlist=None):
   # routine to get output fields from binary files:

   ###########################################################
   # get dimensions and variable names from .b file
   bfile = afile[:-2]+'.b'
   binfo = fn_bfile_info(bfile)

   if vlist is None:
      # get all fields in binary files
      recnos   = binfo['recnos']
   else:
      # check vbls are in binary files
      recnos   = {}
      for vbl in vlist:
         vname = check_names(vbl,binfo['recnos'].keys(),stop=True)
         recnos.update({vbl:binfo['recnos'][vname]})

   nv = binfo['Nrecs']
   nx = binfo['nx']
   ny = binfo['ny']
   if binfo['Norder']==1:
      order = 'fortran'
   else:
      order = 'C'
   ###########################################################

   ###########################################################
   # can now read data from .a file
   import os
   sz       = os.path.getsize(afile)
   fmt_size = int(sz/float(nv*nx*ny)) # 4 for single, 8 for double
   aid      = open(afile,'rb')

   out   = {}
   for vbl in recnos:
      out.update({vbl:get_array(aid,recnos[vbl],nx,ny,order=order,fmt_size=fmt_size)})

   aid.close()
   ###########################################################

   # outputs
   return out,binfo
##############################################################

class grid_info:
   def __init__(self,outdir):
      self.afile = outdir+'/wim_grid.a'
      self.bfile = outdir+'/wim_grid.b'
      for dic in fn_read_general_binary(self.afile):
         # grid_prams, info: both dictionaries
         for key in dic:
            setattr(self,key,dic[key])

      self.dx  = self.scvx.min()
      self.dy  = self.scuy.min()
      return

   def get_vecs_xy(self,units='km',corners=True):

      if units=='km':
         sfac  = 1.e-3
      elif units=='m':
         sfac  = 1.

      if corners:
         x  = self.X[:,0]+.5*self.dx
         y  = self.Y[0,:]+.5*self.dy
         x  = sfac*np.concatenate([[x[0]-self.dx],x])
         y  = sfac*np.concatenate([[y[0]-self.dy],y])
      else:
         x  = sfac*self.X[:,0]
         y  = sfac*self.Y[0,:]

      return x,y

   def plot(self,z=None,zlab=None,pobj=None,zlims=None,colorbar=True,show=True):
      # f  = plt.figure(figsize=[6,6],dpi=50)
      # plt.pcolor(grid_prams.X,grid_prams.Y,ice_fields.icec,cmap=cm.jet,vmax=Vmax,vmin=Vmin)
      from matplotlib import pyplot as plt

      # fontname = 'cursive'
      fontname = 'serif'
      # fontname = 'sans-serif'

      if pobj is None:
         pobj  = plot_object()
      fig   = pobj.fig
      ax    = pobj.ax

      if z is None:
         # plot the land mask
         z     = self.LANDMASK.transpose()
         zlab  = "Land mask"
         zlims = [0.,1.]

      # get corners of grid cells for pcolor
      x,y      = self.get_vecs_xy(units='km',corners=True)

      # check if need to transpose
      nyz,nxz  = z.shape
      if nyz+1==len(x) and nxz+1==len(y):
         z  = z.transpose()

      if zlims is not None:
         vmin  = zlims[0]
         vmax  = zlims[1]
      else:
         vmax  = z.max()
         vmin  = z.min()
         dv    = vmax-vmin

         if dv==0:
            # z is const
            vv    = z.mean()
            vmin  = vmin-.1*vv
            vmax  = vmax+.1*vv
         else:
            vmin  = vmin-.1*dv
            vmax  = vmax+.1*dv
      ################################################


      # make pcolor plot
      P  = ax.pcolor(x,y,z,vmin=vmin,vmax=vmax)
      ax.axes.set_aspect('equal')
      ax.set_xlim((x[0],x[-1]))
      ax.set_ylim((y[0],y[-1]))
      ax.set_xlabel("$x$, km", fontsize=16,fontname=fontname)
      ax.set_ylabel("$y$, km", fontsize=16,fontname=fontname)

      ############################################################
      #colorbar:
      if colorbar:
         from matplotlib import ticker
         # only have colorbar if not constant
         # - doesn't work on linux otherwise
         cbar  = fig.colorbar(P)#,ticks=np.arange(0,1+dc,dc))
         if zlab is not None:
            cbar.set_label(zlab, size=14)
         cbar.ax.tick_params(labelsize=14)
         #plt.locator_params(nbins=4)
         #
         cpos     = cbar.ax.get_position().extents
         cpos[2]  = cpos[2]+.15  # cbar width
         cpos[1]  = cpos[1]+.21  # lower height
         cpos[3]  = cpos[3]*.38  # colorbar height
         cbar.ax.set_position(cpos)         

         tick_locator = ticker.MaxNLocator(nbins=5)
         cbar.locator = tick_locator
         cbar.update_ticks()
      ############################################################

      # fonts of axes:
      for tick in P.axes.xaxis.get_major_ticks():
         tick.label.set_fontsize(14)
         tick.label.set_fontname(fontname)
      for tick in P.axes.yaxis.get_major_ticks():
         tick.label.set_fontsize(14)
         tick.label.set_fontname(fontname)

      if show:
         plt.show(fig)

      return pobj
   #######################################################################
