function plot_steps_mesh(rootdir,variables,plot_options);
%% CALL: plot_steps_mesh(rootdir,variables,plot_options);
%% variables is a cell with strings of variables to plot

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

% default plot options
visible     = 0;
apply_mask  = 1;
save_figure = 0;
if ~exist('plot_options','var');
   plot_options.visible       = visible;
   plot_options.apply_mask    = apply_mask;
   plot_options.save_figure   = save_figure;
end
if ~isfield(plot_options,'visible')
   plot_options.visible = visible;
end
if ~isfield(plot_options,'apply_mask')
   plot_options.apply_mask = apply_mask;
end
if ~isfield(plot_options,'save_figure')
   plot_options.apply_mask = save_figure;
end

outdir   = [rootdir,'/mesh']
figdir   = [rootdir,'/figs']
eval(['!mkdir -p ',figdir]);
figdir   = [figdir,'/mesh']
eval(['!mkdir -p ',figdir]);

%%variables to plot
if exist('variables','var');
   vbls  = variables;
   if ~iscell(vbls)
      error('input ''variables'' should be a cell containing variables to be plotted (as strings)');
   end
else
   vbls  = {'Dfloe' ,...            %1
            'Stresses',...          %2
            'Concentration',...     %3
            'Thickness'  ,...       %4
            'Damage',...            %5
            'Nfloes'...             %6
            'M_wind'...             %7
            };

   if 1
      %% shorten default list of var's
      %jkeep = [3,4,5,7];
      jkeep = [1,2,4,5];
      vbls  = vbls(jkeep);
   end
end
Nv = length(vbls);
%cmaps = {};
%lims  = {[0,350]};

% steps to be loaded
dir0  = dir([outdir,'/mesh_*.dat']);
N0    = length(dir0)-2;%%starts from 0, mesh_1000.dat is final state

region_of_zoom = [];
%region_of_zoom = 'framstrait';
is_sequential  = 1;

for step=0:N0

   for k=1:Nv
      vbl   = vbls{k};
      %cmap  = cmaps{k};
      %lim   = lims{k};

      %% initial/final filenames
      fig_full = [figdir,'/',vbl,'/',vbl,'_',num2str(step,'%2.2d'),'.png'];

      if ~(exist(fig_full)&~OVER_WRITE)
         plot_nextsim_c(vbl,step,region_of_zoom,is_sequential,outdir,plot_options);
         %%
         eval(['!mkdir -p ',figdir,'/',vbl]);
         disp(['saving to ',fig_full]);
         %%
         saveas(gcf,fig_full);
         close;
      end%check if fig is present already
   end%loop over variables
end%loop over time steps
