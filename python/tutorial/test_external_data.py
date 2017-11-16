import nextsim_classes as nsc

ncf   = '/Data/sim/tim/data/SWARP_WW3_ARCTIC-12K_20151210.nc'

if 0:
   nbi   = nsc.file_list(rootdir='out_cpp/mesh/',step1=0,step2=0)
   nbi.plot_external_data(ncf,'hs',time_index=0)
else:
   import os
   msf   = os.getenv('SIMDATADIR')+'/data/mesh/wim_grid_FS_5km_split2.msh'
   mso   = nsc.gmsh_mesh(msf)
   mso.plot_external_data(ncf,'hs')
