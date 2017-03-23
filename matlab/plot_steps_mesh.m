function plot_steps_mesh(rootdir);

if ~exist('rootdir','var');
   %% location of outputs
   run_no   = 2;
   if 0
      %%johansen
      rootdir  = '/Volumes/sim/tim';
   else
      %%external hard disk
      rootdir  = '/Volumes/Tim_Ext_HD2/WORK'
   end
   rootdir  = [rootdir,'/Model-Results/neXtSIM/Oban-test16/run',num2str(run_no)];
end
OVER_WRITE  = 0;

outdir   = [rootdir,'/mesh']
figdir   = [rootdir,'/figs']
eval(['!mkdir -p ',figdir]);
figdir   = [figdir,'/mesh']
eval(['!mkdir -p ',figdir]);

%%variables to plot
vbls  = {'Dfloe' ,...            %1
         'Stresses',...          %2
         'Concentration',...     %3
         'Thickness'  ,...       %4
         'Damage',...            %5
         'Nfloes'...             %6
         'M_wind'...             %7
         };
%cmaps = {};
%lims  = {[0,350]};

% steps to be loaded
dir0  = dir([outdir,'/mesh_*.dat']);
N0    = length(dir0)-2;%%starts from 0, mesh_1000.dat is final state

%% ==================================================
%% shorten list of var's
%% check simul_in to see if waves are present
jkeep = 1:7;
if 1
   %%shorten more manually
   %jkeep = [3,4,5,7];
   jkeep = 6;
end
vbls  = vbls(jkeep);
%cmaps = cmaps(jkeep);
%lims  = lims(jkeep);
Nv = length(vbls);
%% ==================================================

domain         = 'topaz';
region_of_zoom = 'framstrait';
is_sequential  = 1;
for step=0:N0

   for k=1:Nv
      vbl   = vbls{k};
      %cmap  = cmaps{k};
      %lim   = lims{k};

      %% initial/final filenames
      fig_full = [figdir,'/',vbl,'/',vbl,'_',num2str(step,'%2.2d'),'.png'];

      if ~(exist(fig_full)&~OVER_WRITE)
         plot_nextsim_c(vbl,step,domain,region_of_zoom,is_sequential,outdir);
         %%
         eval(['!mkdir -p ',figdir,'/',vbl]);
         disp(['saving to ',fig_full]);
         %%
         saveas(gcf,fig_full);
         close;
      end%check if fig is present already
   end%loop over variables
end%loop over time steps
