function main_script_trajectories(...
    date_start_str,...
    date_end_str,...
    saved_simul_in,...
    test_name...
    )

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

% selection of the maximum and minimum time interval
% data are usualy computed over intervals of more or less 3 days
min_deltat=0;    % minimum time interval in days (default = 1 days)
max_deltat=500;    % maximum time interval in days (default = 6 days)


% definition of the dates
date_start=datenum(date_start_str)
date_end  =datenum(date_end_str)

mod_date_start=date_start;
mod_date_end  =date_end  ;


% model flag and load simul_in
model=0;
if(~isempty(saved_simul_in))
    model=1;
    simul_in=read_simul_in(saved_simul_in); 
    
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
  
    %*************************************************************************%
    %*************************************************************************%
    
    if(model==1)
        % load simul_in
        simul_in=read_simul_in(saved_simul_in); 

        directory_simul_out=fileparts(which(saved_simul_in));
        disp(directory_simul_out);
        
        % steps to be loaded
        output_timestep_days=1./simul_in.output_per_day;
        
        round_dnum_1=floor((date_start-simul_in.time_init)/output_timestep_days)*output_timestep_days;
        round_dnum_2=ceil((date_end  -simul_in.time_init)/output_timestep_days)*output_timestep_days;
        
        step_from = round_dnum_1/output_timestep_days;
        step_to   = round_dnum_2/output_timestep_days;
        
        if(step_to==-1)
            error('No simulated files found')
        end
        
        % compute the avalaible period from the simulation      
        from_date= simul_in.time_init+step_from*output_timestep_days; % date corresponding to from_step of the simulation
        to_date= simul_in.time_init+step_to*output_timestep_days; % date corresponding to from_step of the simulation
        disp(['data from the model loaded from ',datestr(from_date),' to ',datestr(to_date)])
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
mod_step_start=1;
mod_step_end=2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
% Loading model data
%-------------------------------------------------------------------------%


% to load results from X-SIM

% load the step 0
step_from
[mesh_out,data_out] = neXtSIM_bin_revert('',[], step_from);

mean_resol_mod=7;%sqrt(mean(mesh.element.surf))/1000;              % mean resolution, the resolution needs to be uniform

% Limit the domain
%    limit_domain_mod={};
%
%    limit_domain_mod{1}.x=[xmin_mod,xmax_mod,xmax_mod,xmin_mod,xmin_mod]*1000;
%    limit_domain_mod{1}.y=[ymin_mod,ymin_mod,ymax_mod,ymax_mod,ymin_mod]*1000;

%     limit_domain_mod(1).x=[-2000,-250,0,-1000,-2000];
%     limit_domain_mod(1).y=[0,-250,1250,1250,0];

load('RGPS_2007-12-01_2008-06-01_traj.mat')

for i=1:length(out)
    tmp_out=out{i};
    m_proj('Stereographic','lon',-45,'lat',90,'radius',60);
    
    lat=[tmp_out.meta.s_w_lat,tmp_out.meta.s_e_lat,tmp_out.meta.n_e_lat,tmp_out.meta.n_w_lat,tmp_out.meta.s_w_lat];
    lon=[tmp_out.meta.s_w_lon,tmp_out.meta.s_e_lon,tmp_out.meta.n_e_lon,tmp_out.meta.n_w_lon,tmp_out.meta.s_w_lon];
    
    if i==1
        figure(2)
        figure(1)
        worldmap([65 90], [180 360]);
        geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
    end
    figure(1)
    geoshow(lat,lon)
    
    [x,y]=m_ll2xy(lon,lat);
    x_corner_stream=(x*6378.273);
    y_corner_stream=(y*6378.273);
    
    figure(2)
    plot(x_corner_stream,y_corner_stream);
    hold on
    
    limit_domain_mod(i).x=x_corner_stream;
    limit_domain_mod(i).y=y_corner_stream;
end
% add an additional domain to avoid points near the coasts 
limit_domain_mod(i).x=[-2000,-1600,-1000,200,700,700,500,-300,-1400,-2000];
limit_domain_mod(i).y=[0,1100,1400,1200,400,-200,-600,-500,-200,0];

figure(2)
plot(limit_domain_mod(i).x,limit_domain_mod(i).y);
hold on

% Fill in with each step
nb_step=step_to-step_from+1;

% we load only the first and last
step_ids=[1,nb_step];

for i=1:length(step_ids)
    
    % load the step i
    [mesh_out,data_out] = neXtSIM_bin_revert('',[], step_from+step_ids(i)-1);
    
    date_tot(i)=from_date+(step_ids(i)-1)*output_timestep_days;
    
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
    
    %%% nodes
    UM{i}=zeros(2*length(mesh_out.Nodes_x(ind_limit_concentration==1)),1);     
    node_id{i}=mesh_out.id(ind_limit_concentration==1);
    node_x{i}=mesh_out.Nodes_x(ind_limit_concentration==1);
    node_y{i}=mesh_out.Nodes_y(ind_limit_concentration==1);
    
    % Only select te nodes that are within limit_domain_mod
    if(~isempty(limit_domain_mod))
        nb_subdomains=length(limit_domain_mod);
        inside_subdomain=zeros(length(node_x{i}),1);
        
        x_target=node_x{i}/1000;
        y_target=node_y{i}/1000;
        for j=1:nb_subdomains
            sign_cross_product=zeros(length(node_x{i}),1);
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
end

step_from=1;
step_to=2;

% load the node_x, node_y and UM for step_start and step_end
node_x_start=node_x {step_from};
node_y_start=node_y {step_from};
id_start    =node_id{step_from};
UM_start    =UM     {step_from}; % UM is the position relative to the actual mesh

node_x_end  =node_x {step_to  };
node_y_end  =node_y {step_to  };
id_end      =node_id{step_to  };
UM_end      =UM     {step_to  }; % UM is the position relative to the actual mesh

% select the nodes present in the two steps
[C,i_start,i_end] = intersect(id_start,id_end);

% compute the position of the nodes
x_start=node_x_start(i_start)+UM_start(2*i_start-1);
y_start=node_y_start(i_start)+UM_start(2*i_start  );

x_end  =node_x_end  (i_end  )+UM_end  (2*i_end  -1);
y_end  =node_y_end  (i_end  )+UM_end  (2*i_end    );


traj_mod.x_start=x_start;
traj_mod.y_start=y_start;
traj_mod.x_end=x_end;
traj_mod.y_end=y_end;

name_output='traj';
data_type='mod';
save([name_output '_' data_type '_' test_name, '_' num2str(mod_date_start), '_' num2str(mod_date_end) '.mat'],'traj_mod')


end

