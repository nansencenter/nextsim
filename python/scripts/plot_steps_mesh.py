import matplotlib
matplotlib.use('Agg')
import os,sys
from getopt import getopt
from datetime import datetime as DT
import numpy as np
import nextsim_classes as nsc

helpmsg  = '\npython plot_steps_mesh.py [rootdir] [v1] [v2] ... [vN]'

# ==========================================================================
# options
FCtype   = 'ice_only'

# defaults
kwargs   = {'figdir':None,\
            'clim1':None,\
            'clim2':None,\
            'step1':None,\
            'step2':None,\
            'date1':None,\
            'date2':None}

# descriptions
dhelp    = {'figdir':': where to save results',\
            'clim1':': lower range for colorbar',\
            'clim2':': upper range for colorbar',\
            'step1':': step to start from',\
            'step2':': step to finish with',\
            'date1':': date to start with (yyyymmdd)',\
            'date2':': date to finish with (yyyymmdd)'}

# defaults, types
lst_str   = ['figdir'] # strings
lst_int   = ['step1','step2'] # int's
lst_dto   = ['date1','date2'] # dates
lst_dbl   = ['clim1','clim2'] # doubles

kw = []
for key in lst_str:
   kw.append(key+'=')
   helpmsg += ' --'+key+'=...'
for key in lst_int:
   kw.append(key+'=')
   helpmsg += ' --'+key+'=...'
for key in lst_dbl:
   kw.append(key+'=')
   helpmsg += ' --'+key+'=...'
for key in lst_dto:
   kw.append(key+'=')
   helpmsg += ' --'+key+'=...'


helpmsg += '\n\n*rootdir is path to results dir'
helpmsg += '\n*v1, v2,..., vN are variable names'
helpmsg += '\nOptions:'
for key in lst_str:
   helpmsg += '\n*'+key+dhelp[key]
for key in lst_int:
   helpmsg += '\n*'+key+dhelp[key]
for key in lst_dbl:
   helpmsg += '\n*'+key+dhelp[key]
for key in lst_dto:
   helpmsg += '\n*'+key+dhelp[key]

args  = []
for ARG in sys.argv[1:]:
   # print(ARG)
   if '--' in ARG and '=' in ARG:
      opt,arg = ARG.split('=')
      # print(opt,arg)
   else:
      # print('KEEP '+ARG)
      args.append(ARG)
      continue

   # strings
   for key in lst_str:
      if opt=='--'+key:
         if arg!='None':
            kwargs[key] = arg
         continue

   # integers
   for key in lst_int:
      if opt=='--'+key:
         if arg!='None':
            kwargs[key] = int(arg)
         continue

   # doubles
   for key in lst_dbl:
      if opt=='--'+key:
         if arg!='None':
            kwargs[key] = np.double(arg)
         continue

   # dates
   for key in lst_dto:
      if opt=='--'+key:
         if arg!='None':
            kwargs[key] = DT.strptime(arg,'%Y%m%d')
         continue
# ==========================================================================

if len(args)<2:
   helpmsg  = '\n\nNot enough arguments.\n\nUsage:'+helpmsg+'\n'
   raise ValueError(helpmsg)
else:
   rootdir  = args[0]   # path to results
   vlist    = args[1:]  # list of variables to plot

if 0:
   # test
   print(rootdir)
   print(vlist)
   print(kwargs)
   raise ValueError(helpmsg)

# ==========================================================================
# initialise file list object

# move some keywords to plotting function
kwplot   = {}
for kw in ['figdir','clim1','clim2']:
   kwplot.update({kw:kwargs[kw]})
   del(kwargs[kw])
flist = nsc.file_list(rootdir,**kwargs)
# ==========================================================================

# ==========================================================================
# make plots

# change default figdir now we know rootdir
if kwplot['figdir'] is None:
   kwplot['figdir']  = rootdir+'/figs'

# colorbar lims
if (kwplot['clim1'] is not None) and (kwplot['clim2'] is not None):
   kwplot.update({'clim':[kwplot['clim1'],kwplot['clim2']]})
else:
   kwplot.update({'clim':None})
for key in ['clim1','clim2']:
   del(kwplot[key])

for vbl in vlist:
   flist.plot_var_all(vbl,**kwplot)
# ==========================================================================
