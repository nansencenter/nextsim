function main_script(...
    date_start_str,...
    date_end_str,...
    saved_simul_in,...
    obs,...
    compute_deformation)

% Example:
% main_script('07-Mar-2008 00:00:00','16-Mar-2008 00:00:00','simul_in_bigarctic10km_test1',1,0)

% The output are saved in defo_mod,... 
% and can be used for plotting or comparison.
% Example of some plotting functions
% plots_vel_map({'defo_mobs_test1_rgps.mat','defo_rgps.mat'},{'Model','obs'},'arctic',1);
% compare({'defo_mobs_test1_rgps.mat','defo_rgps.mat'},{'Model','obs'})

% All the available observations between date_start and date_end are
% treated
% The modeled displacement is computed from mod_date_start to mod_date_end
% mod_date_start to mod_date_end are also used in regridding to identify 
% the observations that are the nearest in time to (mod_date_start+mod_date_end)/2
% model=1;  % 0=model not used, 1=EB, 2=TOPAZ
% obs=1;    % 0=observation mot used, 1=RGPS, 2=GLOBICE, 3=MEASURE (also named EGPS), 4=SAR (Denis), 5=kwok_rgps, 6=RGPS+GLOBICE+MEASURE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% flags
perform_coarsegrainer=0*compute_deformation;  % The coarse graining is used for the scaling analysis
use_drifters = 0 ;                          % use the drifters instead of the nodes

% selection of the maximum and minimum time interval 
% data are usualy computed over intervals of more or less 3 days
min_deltat=0;    % minimum time interval in days (default = 1 days)                        
max_deltat=500;    % maximum time interval in days (default = 6 days) 


% definition of the dates
date_start=datenum(date_start_str)
date_end  =datenum(date_end_str)

if obs==0
    mod_date_start=date_start;
    mod_date_end  =date_end  ;
else
    default_deltat=3; % in days
    if((date_end-date_start)<default_deltat)
        warning('the choosen perdiod is shorter than the default deltat')
        default_deltat=(date_end-date_start);
        warning(['deltat is reduced to ' num2str(default_deltat)])
    end
    date_middle=(date_start+date_end)/2;
    mod_date_start=date_middle-default_deltat/2;
    mod_date_end=date_middle+default_deltat/2;
end

% model flag and load simul_in
model=0;
if(~isempty(saved_simul_in))
    model=1;
    load(saved_simul_in)
    
    % check if the choosen period is covered by the simulation
    days_in_sec=24*3600;
    simul_in_time_init=simul_in.time_init;
    simul_in_time_end=simul_in.time_init+simul_in.duration/days_in_sec;
    if(simul_in_time_init>min(date_start,mod_date_start))
        disp(datestr(simul_in_time_init))
        warning('The simulation may have started after the beginning of the choosen period')
    end
    if(simul_in_time_end<max(date_start,mod_date_start))
        disp(datestr(simul_in_time_end))
        warning('The simulation may have ended before the end of the choosen period')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of parameters selection
% The code below should not change frequently (add flag above when necessary)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
% Load the model output
%-------------------------------------------------------------------------%

if(model)
    simul_in_domain=simul_in.domain;
    simul_in_resol=simul_in.resol;
    simul_in_name=simul_in.simul_in_name;
    
    meshfile=[simul_in_domain simul_in_resol '.mat'];
%     load(meshfile);
    %*************************************************************************%
    %*************************************************************************%
    
    if(model==1)
        % load simul_in
        saved_simul_in =['simul_in_'  meshfile(1:end-4) '_' simul_in_name '.mat'];
        load(saved_simul_in);
        disp(saved_simul_in);
        
        directory_simul_out=fileparts(which(saved_simul_in));
        disp(directory_simul_out);
        
        % steps to be loaded
        output_timestep_days=simul_in.output_timestep/(24*3600);
                    
        round_dnum_1=floor((date_start-simul_in.time_init)/output_timestep_days)*output_timestep_days;
        round_dnum_2=ceil((date_end  -simul_in.time_init)/output_timestep_days)*output_timestep_days;
        
        step_from = round_dnum_1/output_timestep_days;
        step_to   = round_dnum_2/output_timestep_days;
        
        if(step_to==-1)
            error('No simulated files found')
        end
        
        % compute the avalaible period from the simulation
        output_timestep=simul_in.output_timestep;
        output_timestep_days=output_timestep/(24*3600);
        
        from_date= simul_in.time_init+step_from*output_timestep_days; % date corresponding to from_step of the simulation
        to_date= simul_in.time_init+step_to*output_timestep_days; % date corresponding to from_step of the simulation
        disp(['data from the model loaded from ',datestr(from_date),' to ',datestr(to_date)])
        
    elseif(model==2) % topaz (to be tested)
        %*****************************************************************%
        %*****************************************************************%
        load(topaz_file);
        %*****************************************************************%
        %*****************************************************************%
        
        % information about the simulation
        from_step=0; % first output to be considered for the analysis (usualy 0, which means that we also analyse the spinup period)
        step_to=length(ice.time)-1;  % last output to be considered for the analysis
        
        % compute the avalaible period from the simulation
        output_timestep=24*3600;
        output_timestep_days=output_timestep/(24*3600);
        
        % trick to change the date
        time_init=datenum(ice.time_start);
        
        from_date= time_init+from_step*output_timestep_days; % date corresponding to from_step of the simulation
        to_date= time_init+step_to*output_timestep_days; % date corresponding to from_step of the simulation
    else
        error('model not defined')
    end
else
    output_timestep=0;
    output_timestep_days=0;
end

%-------------------------------------------------------------------------%
% Snapshot selection for the model on its grid and to choose when several
% streams are availables
%-------------------------------------------------------------------------%
if(model)
    mod_step_start=round((mod_date_start-from_date)/output_timestep_days)+1;
    mod_step_end=round((mod_date_end-from_date)/output_timestep_days)+1;
    mod_deltat = (mod_step_end-mod_step_start)*output_timestep_days;
    
    if((mod_date_start<from_date))
        error(['the simulated period does not cover the period choosen for the analysis of the model, ' datestr(mod_date_start) ' is before ' datestr(from_date)])
    end
    if((mod_date_end>to_date))
        error(['the simulated period does not cover the period choosen for the analysis of the model, ' datestr(mod_date_end) ' is after ' datestr(to_date)])
    end
    
    if(mod_deltat<min_deltat || mod_deltat>max_deltat)
        error(['the period choosen for the model on its grid gives a deltat ' num2str(mod_deltat) ' out of [' num2str(min_deltat) ',' num2str(max_deltat) ']' ]);
    end
    
    % A trick to be able to load only the 2 simul_out 
    % when there is no comparison to observations
    if(obs==0)
        mod_step_start=1;
        mod_step_end=2;
    end
else
    mod_step_start=0;
    mod_step_end=0;
    mod_deltat=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
% Loading model data
%-------------------------------------------------------------------------%
x_tot=[];
y_tot=[];
date_tot=[];
node_id=[];
node_x=[];
node_y=[];
UM=[];

if(model)
    if(use_drifters && simul_in.drifters==1)
        % only valid when drifters are deployed with a regular spacing, it does not work when drifters are placed according to IABP buoys)
        method=2; % derive deformation from the drifters 
    else
        method=3; % derive deformation from the nodes
    end
else
    method=0; % method is not used when we do not use model output
end

defo_mod=[];
ebox_mod=[];

defo_mobs=[];
ebox_mobs=[];

if(model==1)
% to load results from X-SIM

    % load the step 0
    saved_simul_out=['simul_out_' meshfile(1:end-4) '_' simul_in_name '_step' num2str(step_from) '.mat'];
    load(saved_simul_out);
  
    % Load/read mesh
    if isfield(simul_out, 'bamg')
        % We're using bamg, so the mesh information is stored in the simul_out file
        [mesh] = importbamg(simul_out.bamg.mesh, simul_out.bamg.geom);
        flag_boundary_fix=1;
        flag_boundary_free=0; 
    else
        load(meshfile);
        flag_boundary_fix=simul_in.flag_fix;
        flag_boundary_free=simul_in.flag_free; 
    end  
    
    mean_resol_mod=sqrt(mean(mesh.element.surf))/1000;              % mean resolution, the resolution needs to be uniform

    xmin_mod=min(mesh.node.x);    % min,max x,y elements coordinates
    ymin_mod=min(mesh.node.y);    % min,max x,y elements coordinates
    xmax_mod=max(mesh.node.x);    % min,max x,y elements coordinates
    ymax_mod=max(mesh.node.y);    % min,max x,y elements coordinates
    
    % Limit the domain 
%    limit_domain_mod={};
% 
%    limit_domain_mod{1}.x=[xmin_mod,xmax_mod,xmax_mod,xmin_mod,xmin_mod]*1000;
%    limit_domain_mod{1}.y=[ymin_mod,ymin_mod,ymax_mod,ymax_mod,ymin_mod]*1000;
  
    limit_domain_mod(1).x=[-2000,-250,0,-1000,-2000];
    limit_domain_mod(1).y=[0,-250,1250,1250,0];

%     load('RGPS_2007-12-01_2008-06-01_traj.mat')
% 
%     for i=1:length(out)
%         tmp_out=out{i};
%         m_proj('Stereographic','lon',-45,'lat',90,'radius',60);
% 
%         lat=[tmp_out.meta.s_w_lat,tmp_out.meta.s_e_lat,tmp_out.meta.n_e_lat,tmp_out.meta.n_w_lat,tmp_out.meta.s_w_lat];
%         lon=[tmp_out.meta.s_w_lon,tmp_out.meta.s_e_lon,tmp_out.meta.n_e_lon,tmp_out.meta.n_w_lon,tmp_out.meta.s_w_lon];
% 
%         if i==1
%             figure(2)
%             figure(1)
%             worldmap([65 90], [180 360]);
%             geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
%         end
%         figure(1)
%         geoshow(lat,lon)
% 
%         [x,y]=m_ll2xy(lon,lat);
%         x_corner_stream=(x*6378.273);
%         y_corner_stream=(y*6378.273);
% 
%         figure(2)
%         plot(x_corner_stream,y_corner_stream);
%         hold on
% 
%         limit_domain_mod(i).x=x_corner_stream;
%         limit_domain_mod(i).y=y_corner_stream;
%     end
    
    

    
    % Fill in with each step
    nb_step=step_to-step_from+1;

    if obs==0,
        % we load only the first and last 
        step_ids=[1,nb_step];
    else
        step_ids=1:nb_step;
    end
    for i=1:length(step_ids)
        
        % load the step i
        saved_simul_out=['simul_out_' meshfile(1:end-4) '_' simul_in_name '_step' num2str(step_from+step_ids(i)-1) '.mat'];
        load(saved_simul_out);
        
        date_tot(i)=from_date+(step_ids(i)-1)*output_timestep_days;
        
        %%% nodes
        UM{i}=simul_out.UM(simul_out.ind_reduced);     
                
        node_reduced=(simul_out.ind_reduced(1:2:end-1)+1)/2;
        node_id{i}=simul_out.node_id(node_reduced);
        node_x{i}=simul_out.bamg.mesh.Vertices(node_reduced,1)*1000;
        node_y{i}=simul_out.bamg.mesh.Vertices(node_reduced,2)*1000;
        
        
        % Only select te nodes that are within limit_domain_mod
        if(~isempty(limit_domain_mod))
            nb_subdomains=length(limit_domain_mod);
            inside_subdomain=zeros(length(node_reduced),1);
            
            x_target=node_x{i}/1000;
            y_target=node_y{i}/1000;
            for j=1:nb_subdomains
                sign_cross_product=zeros(length(node_reduced),1);
                tmp_subdomain=limit_domain_mod(j);
                for k=1:4,
                    x1=tmp_subdomain.x(k);
                    x2=tmp_subdomain.x(k+1);
                    y1=tmp_subdomain.y(k);
                    y2=tmp_subdomain.y(k+1);
                    sign_cross_product=sign_cross_product+sign((x_target-x1).*(y2-y1)-(y_target-y1).*(x2-x1));
                end
                inside_subdomain=inside_subdomain+(abs(sign_cross_product)==4);
            end
            
            tmp=UM{i};
            ind_inside_subdomain=find(inside_subdomain);
            new_UM=zeros(2*length(ind_inside_subdomain),1);
            new_UM(1:2:end-1)=tmp(2*ind_inside_subdomain-1);
            new_UM(2:2:end  )=tmp(2*ind_inside_subdomain  );
            UM{i}=new_UM;
            
            tmp=node_id{i};
            node_id{i}=tmp(inside_subdomain>0);
            tmp=node_x{i};
            node_x{i}=tmp(inside_subdomain>0);
            tmp=node_y{i};
            node_y{i}=tmp(inside_subdomain>0);
        end
        
        %%% drifters
        if(use_drifters)
            % load the step i
            saved_drifters=['drifters_' meshfile(1:end-4) '_' simul_in_name '_step' num2str(step_from+step_ids(i)-1) '.mat'];
            load(saved_drifters);
            size(drifters.x')
            x_tot(:,i) = drifters.x';
            y_tot(:,i) = drifters.y';
        end
        size(node_y{i})
    end

    
elseif(model>1)
    error('model not defined')
end

%-------------------------------------------------------------------------%
% Loading data
%-------------------------------------------------------------------------%
    
% Stream selection
reject_stream='y'; %'b'; % 'y' to do not reject any stream

clear 'defo_obs'
clear 'ebox_rgps'

xmin_obs=inf;    % min,max x,y elements coordinates
ymin_obs=inf;    % min,max x,y elements coordinates
xmax_obs=-inf;
ymax_obs=-inf;

if(model)
   limit_domain=[xmin_mod,xmax_mod,ymin_mod,ymax_mod]*1000;
else
   limit_domain=[];
end

if(obs==1)
    out=RGPS_drift_deformation(date_start_str,date_end_str,compute_deformation,limit_domain);
    obs_data_ind=struct('ind',find((out.stream~=reject_stream).*((out.dnum-out.deltat)>=date_start).*(out.dnum<=date_end).*(out.deltat>=min_deltat).*(out.deltat<=max_deltat)));
    obs_data=out;
    xmin_obs=min(xmin_obs,min(out.xy(:,1))) ;    
    ymin_obs=min(ymin_obs,min(out.xy(:,2))) ;
    xmax_obs=max(xmax_obs,max(out.xy(:,1))) ;
    ymax_obs=max(ymax_obs,max(out.xy(:,2))) ;
    if(compute_deformation)  mean_resol_obs=sqrt(mean(out.area)) ; else mean_resol_obs=sqrt(10^2/2); end
elseif(obs==2)
    out=GLOBICE_drift_deformation(date_start_str,date_end_str,compute_deformation,limit_domain);
    obs_data_ind=struct('ind',find(((out.dnum-out.deltat)>=date_start).*(out.dnum<=date_end).*(out.deltat>=min_deltat).*(out.deltat<=max_deltat)));
    obs_data=out;
    xmin_obs=min(xmin_obs,min(out.xy(:,1))) ;    
    ymin_obs=min(ymin_obs,min(out.xy(:,2))) ;
    xmax_obs=max(xmax_obs,max(out.xy(:,1))) ;
    ymax_obs=max(ymax_obs,max(out.xy(:,2))) ;
    if(compute_deformation)  mean_resol_obs=sqrt(mean(out.area)) ; else mean_resol_obs=sqrt(5^2/2); end
elseif(obs==3)
    out=EGPS_drift_deformation(date_start_str,date_end_str,compute_deformation,limit_domain);
    obs_data_ind=struct('ind',find(((out.dnum-out.deltat)>=date_start).*(out.dnum<=date_end).*(out.deltat>=min_deltat).*(out.deltat<=max_deltat)));
    obs_data=out;
    xmin_obs=min(xmin_obs,min(out.xy(:,1))) ;    
    ymin_obs=min(ymin_obs,min(out.xy(:,2))) ;
    xmax_obs=max(xmax_obs,max(out.xy(:,1))) ;
    ymax_obs=max(ymax_obs,max(out.xy(:,2))) ;
    if(compute_deformation)  mean_resol_obs=sqrt(mean(out.area)) ; else mean_resol_obs=sqrt(10^2/2); end
elseif(obs==4)
    out=SAR_drift_deformation(date_start_str,date_end_str,compute_deformation,limit_domain);
    obs_data_ind=struct('ind',1:length(out.dnum));
    obs_data=out;
    xmin_obs=min(xmin_obs,min(out.xy(:,1))) ;    
    ymin_obs=min(ymin_obs,min(out.xy(:,2))) ;
    xmax_obs=max(xmax_obs,max(out.xy(:,1))) ;
    ymax_obs=max(ymax_obs,max(out.xy(:,2))) ;
    if(compute_deformation)  mean_resol_obs=sqrt(mean(out.area)) ; else mean_resol_obs=3.9; end % The mean_resol_obs has been computed from the output of SAR_drift_deformation
elseif(obs==5)

    if(~compute_deformation)
        error('The kwok_rgps dataset is only treated for deformation. For the drift, please use RGPS')
    end
    define_default_data_path;

    files=dir([indir, 'RGPS_ice_drift/kwok_def/RGPS_*']);
    for i=1:length(files)
        file_date_start=datenum(files(i).name( 6:15));
        file_date_end =datenum(files(i).name(17:26));
        if((file_date_start<=date_start) && (file_date_end>=date_end))
            infile = [indir 'RGPS_ice_drift/kwok_def/' files(i).name];
        end
    end
    if(isempty(infile))
        error('No available data')
    end
    out=load(infile);
    
    [invar]=invariants(out.dudx,out.dudy,out.dvdx,out.dvdy)
    
    min_surf=50;
    max_surf=200;
    max_eps=1;  % This criterion on the total deformation ids discussed in Stern et al. It seems that is was also used for Marsan et al. 2004
    ind=find((invar.eps<=max_eps).*((out.dnum-out.deltat)>=date_start).*(out.dnum<=date_end).*(out.deltat>=min_deltat).*(out.deltat<=max_deltat).*(out.area>=min_surf).*(out.area<=max_surf).*(out.init_area<=max_surf));
    
    ls=sqrt(out.area(ind))';
    x_center=out.x(ind)';
    y_center=out.y(ind)';
    x_patch=[x_center-ls/2,x_center+ls/2,x_center+ls/2,x_center-ls/2];
    y_patch=[y_center-ls/2,y_center-ls/2,y_center+ls/2,y_center+ls/2];
    
    defo_obs.data.dnum=out.dnum(ind)';
    defo_obs.data.uv=zeros(length(ind),2);
    defo_obs.data.xy_tricorner(:,:,1)=x_patch*1e3;
    defo_obs.data.xy_tricorner(:,:,2)=y_patch*1e3;
    defo_obs.data.x=out.x(ind)'*1e3;
    defo_obs.data.y=out.y(ind)'*1e3;
    defo_obs.data.stream=out.stream(ind)';
    defo_obs.data.deltat=out.deltat(ind)';
    defo_obs.data.area=out.area(ind)'*1e6;
    defo_obs.data.scale=sqrt(out.area(ind)')*1e3;
    defo_obs.data.dudx=out.dudx(ind)';
    defo_obs.data.dudy=out.dudy(ind)';
    defo_obs.data.dvdx=out.dvdx(ind)';
    defo_obs.data.dvdy=out.dvdy(ind)';
    defo_obs.data.indices=zeros(3,length(ind));
    
    mean_resol_mod_kwok=10*2;
    min_area_factor=1;
    [ebox_obs]=regridding(order,xmin,ymin,mean_resol_mod_kwok,defo_obs.data.x/1000,defo_obs.data.y/1000,defo_obs.data.stream,defo_obs.data.dnum,defo_obs.data.deltat,(mod_date_start+mod_date_end)/2);
    
    defo_obs.data_bin=[];
    defo_obs.data_bin_mean=[];
    
    obs_data=[];
    obs_data_ind=[];
elseif(obs==6)
    % RGPS
    out=RGPS_drift_deformation(date_start_str,date_end_str,compute_deformation,limit_domain);
    xmin_obs=min(xmin_obs,min(out.xy(:,1))) ;    
    ymin_obs=min(ymin_obs,min(out.xy(:,2))) ;
    xmax_obs=max(xmax_obs,max(out.xy(:,1))) ;
    ymax_obs=max(ymax_obs,max(out.xy(:,2))) ;
    if(compute_deformation)  mean_resol_obs_1=sqrt(mean(out.area)) ; else mean_resol_obs_1=sqrt(10^2/2); end
    rgps_data_ind=find((out.stream~=reject_stream).*((out.dnum-out.deltat)>=date_start).*(out.dnum<=date_end).*(out.deltat>=min_deltat).*(out.deltat<=max_deltat));
    if(isfield(out,'x')) out=rmfield(out,'x'); end
    if(isfield(out,'y')) out=rmfield(out,'y'); end
    rgps_data=out;
    % GLOBICE
    out=GLOBICE_drift_deformation(date_start_str,date_end_str,compute_deformation,limit_domain);
    xmin_obs=min(xmin_obs,min(out.xy(:,1))) ;    
    ymin_obs=min(ymin_obs,min(out.xy(:,2))) ;
    xmax_obs=max(xmax_obs,max(out.xy(:,1))) ;
    ymax_obs=max(ymax_obs,max(out.xy(:,2))) ;
    if(compute_deformation)  mean_resol_obs_2=sqrt(mean(out.area)) ; else mean_resol_obs_2=sqrt(5^2/2); end
    globice_data_ind=find(((out.dnum-out.deltat)>=date_start).*(out.dnum<=date_end).*(out.deltat>=min_deltat).*(out.deltat<=max_deltat));
    globice_data=out;
    % EGPS
    out=EGPS_drift_deformation(date_start_str,date_end_str,compute_deformation,limit_domain);
    xmin_obs=min(xmin_obs,min(out.xy(:,1))) ;    
    ymin_obs=min(ymin_obs,min(out.xy(:,2))) ;
    xmax_obs=max(xmax_obs,max(out.xy(:,1))) ;
    ymax_obs=max(ymax_obs,max(out.xy(:,2))) ;
    if(compute_deformation)  mean_resol_obs_3=sqrt(mean(out.area)) ; else mean_resol_obs_3=sqrt(10^2/2); end
    egps_data_ind=find(((out.dnum-out.deltat)>=date_start).*(out.dnum<=date_end).*(out.deltat>=min_deltat).*(out.deltat<=max_deltat));
    egps_data=out;
    % Merge the 3 datasets
    mean_resol_obs=mean([mean_resol_obs_1,mean_resol_obs_2,mean_resol_obs_3]) ;
    obs_data=[rgps_data,globice_data,egps_data];
    obs_data_ind=[struct('ind',rgps_data_ind),struct('ind',globice_data_ind),struct('ind',egps_data_ind)];
else
    obs_data=[];
    obs_data_ind=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data preparation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(obs~=5)
    
    %-------------------------------------------------------------------------%
    % definition of the parameters used for the regridding
    %-------------------------------------------------------------------------%
    if(model==0)
        mean_resol=mean_resol_obs;
        xmin=xmin_obs;   
        ymin=ymin_obs;
        xmax=xmax_obs;
        ymax=ymax_obs;
    else
        mean_resol=mean_resol_mod;
        xmin=xmin_mod;   
        ymin=ymin_mod;
        xmax=xmax_mod;
        ymax=ymax_mod;
    end
        
    xmin=xmin-mean_resol;    % min,max x,y elements coordinates
    ymin=ymin-mean_resol;    % min,max x,y elements coordinates
    xmax=xmax+mean_resol;    % min,max x,y elements coordinates
    ymax=ymax+mean_resol;    % min,max x,y elements coordinates    
   
    max_domain_size=max(xmax-xmin,ymax-ymin); % max x or y length of the domain
    
    min_box_size=2*mean_resol;          % min size of the boxes used for the regridding and the coarsegraining
    
    min_area_factor=0.5;
    
    order=ceil(log2(max_domain_size/min_box_size)) %
      
    [defo_obs.data,defo_mod.data,defo_mobs.data,ebox_obs,ebox_mod,ebox_mobs] = ...
        prepare_drift_deformation(...
        order,...
        xmin,...
        ymin,...
        xmax,...
        ymax,...
        mean_resol,...
        min_box_size,...
        obs_data,...
        output_timestep,...
        x_tot,...
        y_tot,...
        date_tot,...
        obs_data_ind,...
        mod_step_start,...
        mod_step_end,...
        mod_date_start,...
        mod_date_end,...
        mod_deltat,...
        model,...
        node_id,...
        node_x,...
        node_y,...
        UM,...
        method,...
        compute_deformation);
end

defo_obs.data.indices=[];
if(~isempty(ebox_obs))
    ind_mask=find(ebox_obs.mask);
    for i=ind_mask,
        defo_obs.data.indices  = [defo_obs.data.indices;ebox_obs.no_overlap{i}];   
    end 
else
    defo_obs.data.indices=[]; 
end

defo_mobs.data.indices=[];
if(~isempty(ebox_mobs))
    ind_mask=find(ebox_mobs.mask);
    for i=ind_mask,
        defo_mobs.data.indices  = [defo_mobs.data.indices;ebox_mobs.no_overlap{i}];   
    end 
else
    defo_mobs.data.indices=[]; 
end

defo_mod.data.indices=[];
if(~isempty(ebox_mod))
    ind_mask=find(ebox_mod.mask);
    for i=ind_mask,
        defo_mod.data.indices  = [defo_mod.data.indices;ebox_mod.no_overlap{i}];   
    end 
else
    defo_mod.data.indices=[]; 
end

if(perform_coarsegrainer)
    %-------------------------------------------------------------------------%
    % RGPS Analysis
    %-------------------------------------------------------------------------%
    [defo_obs.data_bin,defo_obs.data_bin_mean] = coarsegrainer(ebox_obs,defo_obs.data,order+1,min_box_size,min_area_factor);
    
    %-------------------------------------------------------------------------%
    % Model as Observation
    %-------------------------------------------------------------------------%
    [defo_mobs.data_bin,defo_mobs.data_bin_mean] = coarsegrainer(ebox_mobs,defo_mobs.data,order+1,min_box_size,min_area_factor);
    
    %-------------------------------------------------------------------------%
    % Model Analysis
    %-------------------------------------------------------------------------%
    [defo_mod.data_bin,defo_mod.data_bin_mean] = coarsegrainer(ebox_mod,defo_mod.data,order+1,min_box_size,min_area_factor);
end

if(compute_deformation)
    name_output='defo';
else
    name_output='drift';
end

if(obs==1)
    save([name_output, '_' num2str(date_start), '_' num2str(date_end) '_rgps.mat'],'defo_obs','ebox_obs')
    if(model)    save([name_output '_mobs_' simul_in_name, '_' num2str(date_start), '_' num2str(date_end) '_rgps.mat'],'defo_mobs','ebox_mobs'); end
elseif(obs==2)
    save([name_output, '_' num2str(date_start), '_' num2str(date_end) '_globice.mat'],'defo_obs','ebox_obs')
    if(model)    save([name_output '_mobs_' simul_in_name, '_' num2str(date_start), '_' num2str(date_end) '_globice.mat'],'defo_mobs','ebox_mobs'); end
elseif(obs==3)
    save([name_output, '_' num2str(date_start), '_' num2str(date_end) '_egps.mat'],'defo_obs','ebox_obs')
    if(model)    save([name_output '_mobs_' simul_in_name, '_' num2str(date_start), '_' num2str(date_end) '_egps.mat'],'defo_mobs','ebox_mobs'); end
elseif(obs==4)
    save([name_output, '_' num2str(date_start), '_' num2str(date_end) '_SAR.mat'],'defo_obs','ebox_obs')
    if(model)    save([name_output '_mobs_' simul_in_name, '_' num2str(date_start), '_' num2str(date_end) '_SAR.mat'],'defo_mobs','ebox_mobs'); end
elseif(obs==5)
    save([name_output, '_' num2str(date_start), '_' num2str(date_end) '_kwok.mat'],'defo_obs','ebox_obs')
elseif(obs==6)
    save([name_output, '_' num2str(date_start), '_' num2str(date_end) '_obs.mat'],'defo_obs','ebox_obs')
    if(model)    save([name_output '_' simul_in_name, '_' num2str(date_start), '_' num2str(date_end) '_mobs.mat'],'defo_mobs','ebox_mobs'); end
end

if(model~=0)
    % Adding displacement
    node_x = mesh.node.x';
    node_y = mesh.node.y';
    
    %Selecting the mesh boundaries
    boundary = mesh.boundary.from_msh;
    %Selecting only closed boundaries
    b_fix = find(flag_boundary_fix==boundary(:,3));
    %Selecting only closed boundaries
    b_free = find(flag_boundary_free==boundary(:,3));
    
    defo_mod.boundaryX  = node_x(boundary(b_fix ,1:2,1));
    defo_mod.boundaryY  = node_y(boundary(b_fix ,1:2,1));
    defo_mod.boundaryXo = node_x(boundary(b_free,1:2,1));
    defo_mod.boundaryYo = node_y(boundary(b_free,1:2,1));
    
    save([name_output '_mod_' simul_in_name, '_' num2str(mod_date_start), '_' num2str(mod_date_end) '.mat'],'defo_mod','ebox_mod')
end

end

