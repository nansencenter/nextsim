function plot_steps_mesh(rootdir,variables,plot_options);
%% CALL: plot_steps_mesh(rootdir,variables,plot_options);
%% INPUTS:
%% *rootdir is the root directory containing directory called mesh (with results on mesh)
%%  figures are saved to rootdir/figs/[variable name]
%% *variables is a cell with strings of variables to plot
%% *plot_options = [default: NB don't need to specify all fields if want to keep default values]
%%        apply_mask: true           % If true, apply ice mask
%%         plot_grid: 0              % If not zero the mesh lines are plotted. If zoomed out only the mesh lines will be visible
%%   plot_coastlines: 1              % When 1 the actual domain boundaries are plotted, closed in light gray and opened in cyan.
%%                                   %  Note though that plotting the coastlines or the grid makes the figure much heavier
%%         plot_date: 0              % 1 if we want to display the date on the figure
%%         font_size: 14             % Sets font size of colorbar and date
%%  background_color: [0.85,.85,.85] % gray-white color. A substitute could be gray [0.5 0.5 0.5]
%%     figure_format: '-png'         % can be pdf, tiff, png or jpeg
%%       pic_quality: '-r300'        % Resolution for eps, pdf, tiff, and png
%%     show_vec_dirn: 0              % If plotting vector magnitude, show the direction as arrows

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
N0    = -1;
step1 = 1e30;
step2 = -1e30;
for j=1:length(dir0)
   f  = dir0(j).name;
   nf = length(f);
   cstep = f(6:nf-4);
   try
      %% eg don't want 'init'
      step  = str2num(cstep);
   catch ME
      continue
   end

   if step<length(dir0)+5
      %% eg don't want final - usually 1000
      step1 = min(step1,step);
      step2 = max(step2,step);
   end
end

disp(['Plotting steps from ',num2str(step1),' to ',num2str(step2),'...']);
disp(' ');

region_of_zoom = [];
%region_of_zoom = 'framstrait';
is_sequential  = 1;

for step=step1:step2

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
