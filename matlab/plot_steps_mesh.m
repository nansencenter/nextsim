function plot_steps_mesh(rootdir,variables,plot_options,step_range);
%% CALL: plot_steps_mesh(rootdir,variables,plot_options,step_range);
%% INPUTS:
%% *rootdir is the root directory containing results on mesh
%%  figures are saved to rootdir/figs/[variable name]
%% *variables is a cell with strings of variables to plot
%% *plot_options = structure eg.
%%    OVER_WRITE: 1           % overwrite figure even if already present, otherwise skip steps that are already plotted
%%    RESPLOT: 0              % if 1, use resplot.m (faster, and doesn't apply restricted range to variables),
%%                            % else use plot_nextsim_c.m (nicer plots)
%%    figdir: [1x99 char]     % where to save figures (default is [rootdir,'/figs'])
%%
%% plot_options can also be used to pass options to plot_nextsim_c.m
%% - default options used are:
%%        apply_mask: true         % If true, apply ice mask
%%         save_figure: 0          % Don't save fig inside plot_nextsim_c
%%         visible: 0              % Don't show figures while plotting
%%
%% *step_range can be used to restrict files to plot:
%% step_range  = []: plot all result files ['field_',num2str(step),'.bin'] where step is numeric
%% step_range  = [step1 step2]: similar to [], but where step1<=step<=step2
%% step_range  = name_filter (a string), to only get files following pattern ['field_',name_filter,'*.bin']

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

% default options that are not for plot_nextsim_c.m
other_prams.OVER_WRITE  = 0;
other_prams.RESPLOT     = 0;
other_prams.figdir      = [rootdir,'/figs'];%%default place to save figures

% default plot options for plot_nextsim_c.m
po_def.visible       = 0;
po_def.apply_mask    = 1;
po_def.save_figure   = 0;

if ~exist('step_range','var'); step_range  = []; end
if ~ischar(step_range)
   %%look for ['mesh_*.dat']
   name_filter = '';
else
   %%look for ['mesh_', name_filter, '*.dat']
   name_filter = step_range;
   step_range  = [];
end

% ========================================================================
if ~exist('plot_options','var');
   plot_options   = po_def;
else
   if isempty(plot_options)
      plot_options   = po_def;
   else
      fields   = fieldnames(po_def);
      for n=1:length(fields)
         fld   = fields{n};
         if ~isfield(plot_options,fld)
            plot_options.(fld)   = po_def.(fld);
         end
      end
   end
end

%% can also pass other_prams (options that don't get passed to plot_nextsim_c.m)
%% with plot_options
fields   = fieldnames(other_prams);
for n=1:length(fields)
   fld   = fields{n};
   if isfield(plot_options,fld)
      other_prams.(fld) = plot_options.(fld);
      plot_options      = rmfield(plot_options,fld);
   end
   eval([fld,' = other_prams.',fld,';']);
end
disp('Options (not for plot_nextsim_c):')
disp(other_prams);
disp('Options (for plot_nextsim_c):')
disp(plot_options);
clear other_prams;
% ========================================================================

%outdir   = [rootdir,'/mesh']
outdir   = rootdir;
eval(['!mkdir -p ',figdir]);
figdir   = [figdir,'/mesh'];
eval(['!mkdir -p ',figdir]);
simul_in = read_simul_in([rootdir,'/nextsim.log'],0);

% ========================================================================
%% default variables to plot
%% - standard variables
vbls_def = {
            'Concentration',...     %1
            'Thickness'  ,...       %2
            'Damage'...             %3
            };

if simul_in.simul.use_wim==1
   %% waves-in-ice vars
   tmp   =  {
             'Stress_waves_ice',...  %4
             'Dfloe' ,...            %5
             'Nfloes',...            %6
            };
   vbls_def(end+1:end+length(tmp))  = tmp;

   if simul_in.nextwim.export_diags_mesh==1
      %% diagnostic waves-in-ice vars
      tmp   =  {
               'Stokes_drift',...      %7
               'Hs',...                %8
               'Tp',...                %9
               'MWD',...               %10
               };
      vbls_def(end+1:end+length(tmp))  = tmp;
   end
   clear tmp;
end

if (simul_in.simul.save_forcing_field==1)&(strcmp(simul_in.setup.atmosphere_type,'constant')==0) 
   vbls_def{end+1}   = 'M_wind'; %11
end

if 0
   %% shorten default list of var's
   jkeep = [1,2,3,4,5,7,8];
   vbls_def  = vbls_def(jkeep);
end
% ========================================================================

if exist('variables','var');
   vbls  = variables;
   if length(vbls)==0%[] or {}
      vbls  = vbls_def;
   elseif ~iscell(vbls)
      error(['input ''variables'' should be a cell containing variables to be plotted (as strings)',...
             ' (use default vars: [] or {})']);
   end
   clear variables;
else
   vbls  = vbls_def;
end
Nv = length(vbls);
vbls
clear vbls_def;
%cmaps = {};
%lims  = {[0,350]};

% steps to be loaded
dir0  = dir([outdir,'/mesh_',name_filter,'*.dat']);
if isempty(name_filter)
   %%only consider numeric steps
   step1 = 1e30;
   step2 = -1e30;
   for j=1:length(dir0)
      f  = dir0(j).name;
      nf = length(f);
      cstep = f(6:nf-4);
      %%
      [step,status]  = str2num(cstep);
      if status==0
         continue
      end
      if step<length(dir0)+5
         %% eg don't want final - usually 1000
         step1 = min(step1,step);
         step2 = max(step2,step);
      end
   end
   if ~isempty(step_range)
      step1 = max(step_range(1),step1);
      step2 = min(step_range(2),step2);
   end

   disp(['Plotting steps from ',num2str(step1),' to ',num2str(step2),'...']);
   steps = {};
   for step=step1:step2
      steps{end+1}   = num2str(step);
   end
   Nfmt  = num2str(length(steps{end}));
   fmt   = ['%',Nfmt,'.',Nfmt,'d'];
else
   steps = {};
   for j=1:length(dir0)
      f              = dir0(j).name;
      nf             = length(f);
      cstep          = f(6:nf-4);
      steps{end+1}   = cstep;
   end
end

Nsteps   = length(steps);
disp(' ');

region_of_zoom = [];
%region_of_zoom = 'framstrait';
is_sequential  = 1;

for nstep=1:Nsteps
   step  = steps{nstep};

   for k=1:Nv
      vbl   = vbls{k};
      %cmap  = cmaps{k};
      %lim   = lims{k};
      po_tmp   = plot_options;
      if strcmp(vbl,'Stokes_drift')|...
            strcmp(vbl,'Hs')|...
            strcmp(vbl,'Tp')|...
            strcmp(vbl,'MWD')
         po_tmp.apply_mask = 0;
      end

      %% initial/final filenames
      if isempty(name_filter)
         %%if numeric
         fig_full = [figdir,'/',vbl,'/',vbl,'_',num2str(str2num(step),fmt),'.png'];
      else
         fig_full = [figdir,'/',vbl,'/',vbl,'_',step,'.png'];
      end

      if ~(exist(fig_full)&~OVER_WRITE)
         disp([vbl,' - ',step,' (',num2str(nstep),'/',num2str(Nsteps),')']);
         if ~RESPLOT
            plot_nextsim_c(vbl,step,region_of_zoom,is_sequential,outdir,po_tmp,simul_in);
         else
            resplot(vbl,step,outdir,po_tmp);
         end
         clear po_tmp;
         %%
         eval(['!mkdir -p ',figdir,'/',vbl]);
         disp(['saving to ',fig_full]);
         disp(['-----------------------------------------------------',...
               '-------------------------------']);
         disp(' ');
         %%
         saveas(gcf,fig_full);
         close all;
      end%check if fig is present already
   end%loop over variables
end%loop over time steps
