import nextsim_funs as nsf
import nextsim_classes as nsc

rdir     = '../model/expt_ideal/expt_05/out_cpp/mesh'
logfile  = rdir+'/nextsim.log'
resfile  = rdir+'/field_12.bin'
nbi      = nsc.nextsim_binary_info(resfile,logfile=logfile)
vbl      = 'Concentration'

print('Plotting the mesh')
nbi.plot_mesh()

print('Plotting %s on the mesh' %(vbl))
nbi.plot_var(vbl)

print('Plotting %s interpolated onto a grid' %(vbl))
nbi.imshow(vbl)
