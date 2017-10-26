function [saved_obs_name,saved_mobs_name,saved_mod_name]=main_script_C(...
    date_start_str,...
    date_end_str,...
    saved_simul_in,...
    obs,...
    compute_deformation,...
    limit_domain_model)

% Example:
% main_script('07-Mar-2008 00:00:00','16-Mar-2008 00:00:00','simul_in_bigarctic10km_test1',1,0,0)
% main_script('15-Sept-2007 00:00:00','15-Sept-2007 12:00:00','nextsim.log',0,0,0)

% The output are saved in defo_mod,... 
% and can be used for plotting or comparison.
% Example of some plotting functions
% plots_vel_map({'defo_mobs_test1_rgps.mat','defo_rgps.mat'},{'Model','obs'},'arctic',1,'');
% compare({'defo_mobs_test1_rgps.mat','defo_rgps.mat'},{'Model','obs'})

% All the available observations between date_start and date_end are
% treated
% The modeled displacement is computed from mod_date_start to mod_date_end
% mod_date_start to mod_date_end are also used in regridding to identify 
% the observations that are the nearest in time to (mod_date_start+mod_date_end)/2
% model=1;  % 0=model not used, 1=EB, 2=TOPAZ, 3=neXtSIM
% obs=1;    % 0=observation mot used, 
            % 1=RGPS, 
            % 2=GLOBICE, 
            % 3=MEASURE (also named EGPS), 
            % 4=SAR (Denis), 
            % 5=kwok_rgps, 
            % 6=RGPS+GLOBICE+MEASURE, 
            % 7=DTU sea ice drift
% drift
% limit_domain_model    % 0  take the whole domain, 
%                       % 1 limit the domain to a polygon
%                       % 2 limit the domain to the area covered by the RGPS streams for the winter 2007-2008
%                       % 3 limit the domain to the area covered by the
%                       RGPS streams for the winter 2007-2008 and remove
%                       coastal areas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


RGPS_file_traj='RGPS_2006-12-03_2007-06-01_traj.mat';
starting_date_str=RGPS_file_traj(6:15);

test_name='_test';
test_name='';
%test_name='_viscous';

model_file_traj=['RGPS_Drifters_' starting_date_str(1:4) starting_date_str(6:7) starting_date_str(9:10)];
model_file_traj=[model_file_traj test_name '.nc']



% flags
perform_coarsegrainer=compute_deformation;  % The coarse graining is used for the scaling analysis
use_drifters = 1 ;                          % use the drifters instead of the nodes

% selection of the maximum and minimum time interval 
% data are usualy computed over intervals of more or less 3 days
if(obs==7)
    min_deltat=0.1;  % minimum time interval in days (default = 1 days)                        
    max_deltat=6;    % maximum time interval in days (default = 6 days) 
elseif(obs)
    min_deltat=2.9;    % minimum time interval in days (default = 1 days)                        
    max_deltat=3.1;    % maximum time interval in days (default = 6 days)
else
    min_deltat=0;    % minimum time interval in days (default = 1 days)                        
    max_deltat=1e10; % maximum time interval in days (default = 6 days) 
end   
    
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
    if(strcmp(saved_simul_in(1:8),'simul_in'))
        model=1;
        load(saved_simul_in)
    else
        model=3;
        simul_in=read_simul_in(saved_simul_in);
        simul_in=simul_in.simul;
    end
    
    % check if the choosen period is covered by the simulation
    days_in_sec=24*3600;
    simul_in_time_init=simul_in.time_init;
    disp(simul_in_time_init);
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

if(compute_deformation)
    name_output='defo';
else
    name_output='drift';
end

if(model)
    simul_in_name='test';
end

%-------------------------------------------------------------------------%
% Check if the defo files already exists, if they do, we do nothing.
%-------------------------------------------------------------------------%

saved_obs_name='';
saved_mobs_name='';
saved_mod_name=''; 

if(model)    
    saved_mod_name=[name_output '_mod_' simul_in_name, '_' num2str(mod_date_start), '_' num2str(mod_date_end) '.mat'];
end
if(obs==1)
    saved_obs_name = [name_output, '_' num2str(date_start), '_' num2str(date_end) '_rgps.mat'];
    if(model)
        saved_mobs_name = [name_output '_mobs_' simul_in_name, '_' num2str(date_start), '_' num2str(date_end) '_rgps.mat'];
    end
elseif(obs==2)
    saved_obs_name = [name_output, '_' num2str(date_start), '_' num2str(date_end) '_globice.mat'];
    if(model)    
        saved_mobs_name = [name_output '_mobs_' simul_in_name, '_' num2str(date_start), '_' num2str(date_end) '_globice.mat'];
    end
elseif(obs==3)
    saved_obs_name = [name_output, '_' num2str(date_start), '_' num2str(date_end) '_egps.mat'];
    if(model)    
        saved_mobs_name = [name_output '_mobs_' simul_in_name, '_' num2str(date_start), '_' num2str(date_end) '_egps.mat'];
    end
elseif(obs==4)
    saved_obs_name = [name_output, '_' num2str(date_start), '_' num2str(date_end) '_SAR.mat'];
    if(model)    
        saved_mobs_name = [name_output '_mobs_' simul_in_name, '_' num2str(date_start), '_' num2str(date_end) '_SAR.mat'];
    end
elseif(obs==5)
    saved_obs_name = [name_output, '_' num2str(date_start), '_' num2str(date_end) '_kwok.mat'];
elseif(obs==6)
    saved_obs_name = [name_output, '_' num2str(date_start), '_' num2str(date_end) '_obs.mat'];
    if(model)    
        saved_mobs_name = [name_output '_' simul_in_name, '_' num2str(date_start), '_' num2str(date_end) '_mobs.mat'];
    end
elseif(obs==7)
    saved_obs_name = [name_output, '_' num2str(date_start), '_' num2str(date_end) '_dtu.mat'];
    if(model)    
        saved_mobs_name = [name_output '_mobs_' simul_in_name, '_' num2str(date_start), '_' num2str(date_end) '_dtu.mat'];
    end
end

((exist(saved_mod_name)==2).*(exist(saved_obs_name)==2)*(exist(saved_mobs_name)==2))==1

if(((exist(saved_mod_name)==2).*(exist(saved_obs_name)==2)*(exist(saved_mobs_name)==2))==1)
    return
end



%-------------------------------------------------------------------------%
% Load the model output
%-------------------------------------------------------------------------%

if(model)
    
    % EB in matlab
    if(model==1 || model==3)
        if(model==1)
            simul_in_domain=simul_in.domain;
            simul_in_resol=simul_in.resol;
            simul_in_name=simul_in.simul_in_name;
    
            meshfile=[simul_in_domain simul_in_resol '.mat'];
        end
%         % already loaded
%         % load simul_in
%         saved_simul_in =['simul_in_'  meshfile(1:end-4) '_' simul_in_name '.mat'];
%         load(saved_simul_in);
%         disp(saved_simul_in);
        
        directory_simul_out=fileparts(which(saved_simul_in));
        disp(directory_simul_out);
        
        % steps to be loaded
        if(model==1)
            output_timestep=simul_in.output_timestep;
            output_timestep_days=output_timestep/days_in_sec;
        elseif(model==3)
            output_timestep_days=1/simul_in.output_per_day;
            output_timestep=output_timestep_days*days_in_sec;
            
            if(use_drifters)
                output_timestep_days=simul_in.rgps_drifters_output_time_step;
                output_timestep=output_timestep_days*days_in_sec;
            end
        end
        
        if(use_drifters)
            time_init=datenum(starting_date_str);
        else
            time_init=simul_in.time_init;
        end
            
        if(output_timestep_days~=1)
            round_dnum_1=floor((date_start-time_init)/output_timestep_days)*output_timestep_days;
            round_dnum_2=ceil((date_end  -time_init)/output_timestep_days)*output_timestep_days;
        else
            round_dnum_1=(date_start-time_init);
            round_dnum_2=(date_end  -time_init);    
        end
        date_start
        time_init
        output_timestep_days
        step_from = round_dnum_1/output_timestep_days
        step_to   = round_dnum_2/output_timestep_days
        
        if(step_to==-1)
            error('No simulated files found')
        end
        
        % compute the avalaible period from the simulation   
        from_date= time_init+step_from*output_timestep_days; % date corresponding to from_step of the simulation
        to_date= time_init+step_to*output_timestep_days; % date corresponding to from_step of the simulation
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
    if(use_drifters )
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

limit_domain = define_limit_domain( limit_domain_model );

%-------------- Load model output -------------------
if(model==1)
% to load results from EB ran with the matlab code

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
        
        [mesh, ind_node_fix_bnd, ind_node_free_bnd] = importbamg(simul_out.bamg.mesh, simul_out.bamg.geom);
        
        
        % Filter as a function of the minimum concentration in the
        % surrounding elements
        limit_concentration=0.8;
        node_reduced=(simul_out.ind_reduced(1:2:end-1)+1)/2;
        [C,IA,IB] = intersect(node_reduced,ind_node_free_bnd);
        ind_limit_concentration=ones(length(node_reduced),1);
        nodal_element_connectivity=simul_out.bamg.mesh.NodalElementConnectivity;
        
        for j=node_reduced',
            local_connectivity=nodal_element_connectivity(j,:);
            local_connectivity=local_connectivity(~isnan(local_connectivity));
            min_concentration=min(simul_out.c( local_connectivity));
            if(min_concentration<=limit_concentration)
                ind_limit_concentration(j)=0;
            end
        end
        
        date_tot(i)=from_date+(step_ids(i)-1)*output_timestep_days;
        
        %%% nodes
        selected_nodes_conc=find(ind_limit_concentration==1);
        selected_nodes=zeros(2*length(selected_nodes_conc),1);
        selected_nodes(1:2:end-1)=2*node_reduced(selected_nodes_conc)-1;
        selected_nodes(2:2:end)  =2*node_reduced(selected_nodes_conc)  ;
        
        UM{i}=simul_out.UM(selected_nodes);     
                
        node_reduced=node_reduced(selected_nodes_conc);
        node_id{i}=simul_out.node_id(node_reduced);
        node_x{i}=simul_out.bamg.mesh.Vertices(node_reduced,1)*1000;
        node_y{i}=simul_out.bamg.mesh.Vertices(node_reduced,2)*1000;
        
        %%% drifters
        if(use_drifters)
            % load the step i
            saved_drifters=['drifters_' meshfile(1:end-4) '_' simul_in_name '_step' num2str(step_from+step_ids(i)-1) '.mat'];
            load(saved_drifters);
            size(drifters.x')
            x_tot(:,i) = drifters.x';
            y_tot(:,i) = drifters.y';
        end
        
        
        % Only select te nodes that are within limit_domain
        if(~isempty(limit_domain))
            x_target=node_x{i}/1000;
            y_target=node_y{i}/1000;
            
            inside_subdomain = is_inside_subdomain( limit_domain, x_target, y_target );
            
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
        
    end

elseif(model==3)
% to load results from neXtSIM ran with the C code

    % load the step 0
    step_from
    [mesh_out,data_out] = neXtSIM_bin_revert('',[], step_from);
    
    mean_resol_mod=7;%sqrt(mean(mesh.element.surf))/1000;              % mean resolution, the resolution needs to be uniform
    
    xmin_mod=min(mesh_out.Nodes_x/1000);    % min,max x,y elements coordinates
    ymin_mod=min(mesh_out.Nodes_y/1000);    % min,max x,y elements coordinates
    xmax_mod=max(mesh_out.Nodes_x/1000);    % min,max x,y elements coordinates
    ymax_mod=max(mesh_out.Nodes_y/1000);    % min,max x,y elements coordinates
    
    if(use_drifters) 
        % Load the trajecrtories from the model
        latitude=ncread(model_file_traj,'latitude'); % m
        longitude=ncread(model_file_traj,'longitude'); % m
        time=ncread(model_file_traj,'time'); % m
        time=datenum('1900-01-01 00:00:00')+double(time);
        
        [nb_fake_dimension,nb_buoys_mod,n_obs_mod]=size(latitude);
        
        for i=1:nb_buoys_mod
            latitude(1,i,(latitude(1,i,:)==latitude(1,i,end)).*(longitude(1,i,:)==longitude(1,i,end))>0)=NaN;
            longitude(1,i,(latitude(1,i,:)==latitude(1,i,end)).*(longitude(1,i,:)==longitude(1,i,end))>0)=NaN;
        end
        
        % projection used when reading the RGPS Lagrangian data
        m_proj('Stereographic','lon',-45,'lat',90,'radius',60);
        
        [x_mod,y_mod]=m_ll2xy(squeeze(longitude(1,:,:)),squeeze(latitude(1,:,:)));
        x_mod=double(x_mod)*6378.273;
        y_mod=double(y_mod)*6378.273;
    end
    
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
        [mesh_out,data_out] = neXtSIM_bin_revert('',[], step_from+step_ids(i)-1);
        
        % Filter as a function of the maximum concentration in the
        % surrounding elements
        limit_concentration=0.8;
        ind_limit_concentration=zeros(size(mesh_out.Nodes_x));
        ind_non_empty_element=find(data_out.Concentration>=limit_concentration);
        
        Ne=length(mesh_out.Elements)/3;

        Element=reshape(mesh_out.Elements,[3,Ne]);
        ind_limit_concentration(Element(1,ind_non_empty_element))=1;
        ind_limit_concentration(Element(2,ind_non_empty_element))=1;
        ind_limit_concentration(Element(3,ind_non_empty_element))=1;
        
        date_tot(i)=from_date+(step_ids(i)-1)*output_timestep_days;

        %%% nodes
%         selected_nodes_conc=find(ind_limit_concentration==1);
%         selected_nodes=zeros(2*length(selected_nodes_conc),1);
%         selected_nodes(1:2:end-1)=2*node_reduced(selected_nodes_conc)-1;
%         selected_nodes(2:2:end)  =2*node_reduced(selected_nodes_conc)  ;


        UM{i}=zeros(2*length(mesh_out.Nodes_x(ind_limit_concentration==1)),1);     
        node_id{i}=mesh_out.id(ind_limit_concentration==1);
        node_x{i}=mesh_out.Nodes_x(ind_limit_concentration==1);
        node_y{i}=mesh_out.Nodes_y(ind_limit_concentration==1);
        
        %%% drifters
        if(use_drifters)
            % load the step i
            x_tot(:,i) = x_mod(:,step_from+step_ids(i));
            y_tot(:,i) = y_mod(:,step_from+step_ids(i));
        end
        
        % Only select te nodes that are within limit_domain
        if(~isempty(limit_domain))
            nb_subdomains=length(limit_domain);
            inside_subdomain=zeros(length(node_x{i}),1);
            
            x_target=node_x{i}/1000;
            y_target=node_y{i}/1000;
            for j=1:nb_subdomains
                sign_cross_product=zeros(length(node_x{i}),1);
                tmp_subdomain=limit_domain(j);
                nb_edges_polygon=length(tmp_subdomain.x)-1;
                for k=1:nb_edges_polygon,
                    x1=tmp_subdomain.x(k);
                    x2=tmp_subdomain.x(k+1);
                    y1=tmp_subdomain.y(k);
                    y2=tmp_subdomain.y(k+1);
                    sign_cross_product=sign_cross_product+sign((x_target-x1).*(y2-y1)-(y_target-y1).*(x2-x1));
                end
                inside_subdomain=inside_subdomain+(abs(sign_cross_product)==nb_edges_polygon);
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
        
    end
    
%else
%    error('model not defined')
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
   limit_domain=limit_domain;%[xmin_mod,xmax_mod,ymin_mod,ymax_mod]*1000;
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
elseif(obs==7)
    out=DTU_drift_deformation(date_start_str,date_end_str,compute_deformation,limit_domain);
    obs_data_ind=struct('ind',find(((out.dnum-out.deltat)>=date_start).*(out.dnum<=date_end).*(out.deltat>=min_deltat).*(out.deltat<=max_deltat)));
    obs_data=out;
    xmin_obs=min(xmin_obs,min(out.xy(:,1))) ;    
    ymin_obs=min(ymin_obs,min(out.xy(:,2))) ;
    xmax_obs=max(xmax_obs,max(out.xy(:,1))) ;
    ymax_obs=max(ymax_obs,max(out.xy(:,2))) ;
    if(compute_deformation)  mean_resol_obs=sqrt(mean(out.area)) ; else mean_resol_obs=sqrt(10^2/2); end
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
    elseif(obs==0)
        mean_resol=mean_resol_mod
        xmin=xmin_mod;   
        ymin=ymin_mod;
        xmax=xmax_mod;
        ymax=ymax_mod;
    else
        mean_resol=max(mean_resol_mod,mean_resol_obs)
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

if(obs && model)
    ebox_obs.mask=ebox_obs.mask.*ebox_mobs.mask;
    ebox_mobs.mask=ebox_obs.mask.*ebox_mobs.mask;
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

if(obs==1)
    save(saved_obs_name,'defo_obs','ebox_obs')
    if(model)
        save(saved_mobs_name,'defo_mobs','ebox_mobs'); 
    end
elseif(obs==2)
    save(saved_obs_name,'defo_obs','ebox_obs')
    if(model)    
        save(saved_mobs_name,'defo_mobs','ebox_mobs'); 
    end
elseif(obs==3)
    save(saved_obs_name,'defo_obs','ebox_obs')
    if(model)    
        save(saved_mobs_name,'defo_mobs','ebox_mobs'); 
    end
elseif(obs==4)
    save(saved_obs_name,'defo_obs','ebox_obs')
    if(model)    
        save(saved_mobs_name,'defo_mobs','ebox_mobs'); 
    end
elseif(obs==5)
    save(saved_obs_name,'defo_obs','ebox_obs')
elseif(obs==6)
    save(saved_obs_name,'defo_obs','ebox_obs')
    if(model)    
        save(saved_mobs_name,'defo_mobs','ebox_mobs'); 
    end
elseif(obs==7)
    save(saved_obs_name,'defo_obs','ebox_obs')
    if(model)    
        save(saved_mobs_name,'defo_mobs','ebox_mobs'); 
    end
end

if(model~=0)
    if(model==1)
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
    end
    
    
    save(saved_mod_name,'defo_mod','ebox_mod')
end

end


