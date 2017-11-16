function [obs_data,mod_data,mobs_data,ebox_obs,ebox_mod,ebox_mobs] = ...
    prepare_drift_deformation(...
    order,...
    xmin,...
    ymin,...
    xmax,...
    ymax,...
    mean_resol,...
    min_box_size,...
    out,...
    output_timestep,...
    x_tot,...
    y_tot,...
    date_tot,...
    obs_ind,...
    mod_step_start,...
    mod_step_end,...
    mod_date_start,...
    mod_date_end,...
    mod_deltat,...
    mod,...
    node_id,...
    node_x,...
    node_y,...
    UM,...
    method,...
    compute_deformation)


%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

obs_data = [];
mod_data = [];
mobs_data = [];
ebox_obs = [];
ebox_mod = [];
ebox_mobs = [];

min_size_stream=1;

date_target=(mod_date_start+mod_date_end)/2;

%-------------------------------------------------------------------------%
% RGPS deformation data preparation
%-------------------------------------------------------------------------%
if(~isempty(out))
    
    obs_data.dnum=[];
    obs_data.uv = [];
    obs_data.x = [];
    obs_data.y = [];
    obs_data.stream=[];
    obs_data.deltat=[];
    
    if(compute_deformation)
        obs_data.xy_tricorner=[];
        obs_data.uv_tricorner=[];
        obs_data.area=[];
        obs_data.dudx=[];
        obs_data.dudy=[];
        obs_data.dvdx=[];
        obs_data.dvdy=[];
    end
    
    for k=1:length(out)
               
        obs_data.dnum=[obs_data.dnum;out(k).dnum(obs_ind(k).ind)];
        obs_data.uv = [obs_data.uv;out(k).uv(obs_ind(k).ind,:)*1000]; % m/day
        obs_data.stream=[obs_data.stream;out(k).stream(obs_ind(k).ind)];
        obs_data.deltat=[obs_data.deltat;out(k).deltat(obs_ind(k).ind)];
        
        if(compute_deformation)
            obs_data.xy_tricorner=[obs_data.xy_tricorner;out(k).xy_tricorner(obs_ind(k).ind,:,:)*1000]; % m
            obs_data.uv_tricorner=[obs_data.uv_tricorner;out(k).uv_tricorner(obs_ind(k).ind,:,:)*1000]; % m/day
            obs_data.area=[obs_data.area;out(k).area(obs_ind(k).ind)*1e6]; % m^2
            obs_data.dudx=[obs_data.dudx;out(k).dudx(obs_ind(k).ind)];
            obs_data.dudy=[obs_data.dudy;out(k).dudy(obs_ind(k).ind)];
            obs_data.dvdx=[obs_data.dvdx;out(k).dvdx(obs_ind(k).ind)];
            obs_data.dvdy=[obs_data.dvdy;out(k).dvdy(obs_ind(k).ind)];
        else
            obs_data.x = [obs_data.x;out(k).xy(obs_ind(k).ind,1)*1000];
            obs_data.y = [obs_data.y;out(k).xy(obs_ind(k).ind,2)*1000];
        end
    end
    
    if(compute_deformation)
        obs_data.x = mean(obs_data.xy_tricorner(:,:,1),2);
        obs_data.y = mean(obs_data.xy_tricorner(:,:,2),2);
        obs_data.scale=sqrt(obs_data.area);
    end
    
    [ebox_obs]=regridding(order,xmin,ymin,min_box_size,obs_data.x/1000,obs_data.y/1000,obs_data.stream,obs_data.dnum-obs_data.deltat/2,obs_data.deltat,(mod_date_start+mod_date_end)/2);
end

%-------------------------------------------------------------------------%
% model deformation data preparation
%-------------------------------------------------------------------------%
if(mod)
    
    step_from = mod_step_start;
    step_to   = mod_step_end;
    date_to   = mod_date_end;
    deltat    = mod_deltat;
    
    switch method
        case 2 % derive deformation from the drifters (only valid when drifters are activated)
            
            x_drifters=x_tot(:,step_to)*1000; %m
            y_drifters=y_tot(:,step_to)*1000; %m
            
            u_drifters=(x_tot(:,step_to)-x_tot(:,step_from))/(deltat)*1000; %m/day
            v_drifters=(y_tot(:,step_to)-y_tot(:,step_from))/(deltat)*1000; %m/day
            
            selected_drifters=find(~isnan(x_drifters).*~isnan(y_drifters).*~isnan(u_drifters).*~isnan(v_drifters));
            
            x_tmp=x_drifters(selected_drifters);
            y_tmp=y_drifters(selected_drifters);
            
            u_tmp=u_drifters(selected_drifters); %m/day
            v_tmp=v_drifters(selected_drifters); %m/day
           
        case 3 % derive deformation from the nodes
            
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
            
            % compute the velocity
            u_node=(x_end-x_start)/(deltat); %m/day
            v_node=(y_end-y_start)/(deltat); %m/day
            
            x_tmp=x_end; %m
            y_tmp=y_end; %m
            
            u_tmp=u_node; %m/day
            v_tmp=v_node; %m/day
            
        otherwise
            error('Method not implemented')
    end
    
    % prepare data for def_from_drift
    txy2=[x_tmp,y_tmp]; %m
    trux2=u_tmp; %m/s
    trvy2=v_tmp; %m/s      
    
    % clone of the obs:
    nb_selected_point=length(trux2)
    tdeltat2 =deltat;
    tdnum2   =date_to;
    stream(1:nb_selected_point)='z';
    
    if(compute_deformation)
        
        % parameters to reject bad triangles after the triangulation
        min_area=0;
        max_area=(min_box_size*2)^2*1e6;
        min_minang=0;
        max_long_side=min_box_size*4*1e3;
        showplot=0;
        filter_noise=0;
        min_def=0.02;
        max_level=3;
        
        [ntri,ttdnum,ttdeltat,area,xy_tricorner,uv,uv_tricorner,dudx,dudy,dvdx,dvdy,id_stream,quality_index]=def_from_drift(txy2,tdeltat2,tdnum2,trux2,trvy2,'z',min_area,max_area,min_minang,max_long_side,showplot,filter_noise,min_def,max_level);
        
        if(ntri>0)
            mod_data.dnum=ttdnum*ones(size(area));
            mod_data.uv = uv; % in m/day
            mod_data.xy_tricorner=xy_tricorner; % in m
            mod_data.uv_tricorner=uv_tricorner; % in m/day
            mod_data.stream=id_stream;
            mod_data.deltat=ttdeltat*ones(size(area));
            mod_data.area=area; % km^2
            mod_data.dudx=dudx; % in /day
            mod_data.dudy=dudy; % in /day
            mod_data.dvdx=dvdx; % in /day
            mod_data.dvdy=dvdy; % in /day
        else
            error('empty triangulation')
        end
        
        mod_data.x = mean(mod_data.xy_tricorner(:,:,1),2);
        mod_data.y = mean(mod_data.xy_tricorner(:,:,2),2);
        mod_data.scale=sqrt(mod_data.area);
    else
        mod_data.x = txy2(:,1);
        mod_data.y = txy2(:,2);
        mod_data.dnum=tdnum2*ones(size(mod_data.y));
        mod_data.uv = [trux2,trvy2]; % in m/day
        mod_data.stream='z'*ones(size(mod_data.y));
        mod_data.deltat=tdeltat2*ones(size(mod_data.y));
    end

    [ebox_mod]=regridding(order,xmin,ymin,min_box_size,mod_data.x/1000,mod_data.y/1000);
end

%-------------------------------------------------------------------------%
% Model as if it was from data (mobs)
%-------------------------------------------------------------------------%
if(~isempty(out) && mod)
    
    res=mean_resol/2; % in km
    min_ind_x=floor((xmin-10*res)/res);
    max_ind_x=ceil ((xmax+10*res)/res);
    min_ind_y=floor((ymin-10*res)/res);
    max_ind_y=ceil ((ymax+10*res)/res);
      
    [xq,yq] = meshgrid(min_ind_x*res:res:max_ind_x*res,min_ind_y*res:res:max_ind_y*res); % in km
    
    if(isempty(xq))
        return;
    end
    
    date_target=(mod_date_start+mod_date_end)/2;
    
    mobs_data.dnum=[];
    mobs_data.uv = [];
    mobs_data.stream=[];
    mobs_data.deltat=[];
    mobs_data.x=[];
    mobs_data.y=[];
    mobs_data.scale=[];
    
    if(compute_deformation)
        mobs_data.xy_tricorner=[];
        mobs_data.uv_tricorner=[];
        mobs_data.area=[];
        mobs_data.dudx=[];
        mobs_data.dudy=[];
        mobs_data.dvdx=[];
        mobs_data.dvdy=[];
    end
        
    % for each dataset
    for k=1:length(out)
        
        out_stream=out(k).stream(obs_ind(k).ind);
        out_dnum=out(k).dnum(obs_ind(k).ind);
        out_deltat=out(k).deltat(obs_ind(k).ind);
        
        if(compute_deformation)
            out_xy_tricorner=out(k).xy_tricorner(obs_ind(k).ind,:,:);
        else
            out_xy=out(k).xy(obs_ind(k).ind,:);
        end
        
        % for each stream
        streams=unique(out_stream);
        for s=1:length(streams)
            ind_k=find(out_stream==streams(s));
            
            % for each date
            dnums=unique(out_dnum(ind_k));
            for l=1:length(dnums)
                ind_l=find(out_dnum(ind_k)==dnums(l));
                
                % for each deltat
                deltats=unique(out_deltat(ind_k(ind_l)));
                for m=1:length(deltats)
                    ind_m=find(out_deltat(ind_k(ind_l))==deltats(m));
                    
                    %-------------------------------------------------------------------------%
                    % Relocation
                    % for each rgps point, search its location in the regular grid (? la "fast and furious")
                    %-------------------------------------------------------------------------%
                    
                    % Matching in time
                    output_timestep_days=output_timestep/(24*3600);
                    
                    round_dnum_2=round((dnums(l)-date_tot(1))/output_timestep_days)*output_timestep_days;
                    round_dnum_1=round((dnums(l)-deltats(m)-date_tot(1))/output_timestep_days)*output_timestep_days;
                     
                    step_from = round_dnum_1/output_timestep_days+1;
                    step_to   = round_dnum_2/output_timestep_days+1;
                    date_to   = round_dnum_2; % date_to is not used, we used dnums(1) instead to have exactly the same date as the RGPS
                    deltat=round_dnum_2-round_dnum_1; % deltat is only used to computed the velocity
                                      
                    if(step_from<1)
                        error('The model has not produced results for the choosen period, trick or change the period to fit the simulation')
                    end
                    
                    % compute x_tmp, y_tmp, u_tmp, v_tmp 
                    % either from the drifters or the nodes
                    switch method
                        case 2 % derive deformation from the drifters (only valid when drifters are activated)
                            
                            x_drifters=x_tot(:,step_to)*1000; %m
                            y_drifters=y_tot(:,step_to)*1000; %m
                            
                            u_drifters=(x_tot(:,step_to)-x_tot(:,step_from))/(deltat)*1000; %m/day
                            v_drifters=(y_tot(:,step_to)-y_tot(:,step_from))/(deltat)*1000; %m/day
                            
                            selected_drifters=find(~isnan(x_drifters).*~isnan(y_drifters).*~isnan(u_drifters).*~isnan(v_drifters));
                            
                            x_tmp=x_drifters(selected_drifters);
                            y_tmp=y_drifters(selected_drifters);
                            
                            u_tmp=u_drifters(selected_drifters); %m/day
                            v_tmp=v_drifters(selected_drifters); %m/day
                            
                        case 3 % derive deformation from the nodes
                            
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
                            
                            % compute the velocity
                            u_node=(x_end-x_start)/(deltat); %m/day
                            v_node=(y_end-y_start)/(deltat); %m/day
                            
                            x_tmp=x_end;
                            y_tmp=y_end;
                            
                            u_tmp=u_node; %m/day
                            v_tmp=v_node; %m/day
                            
                        otherwise
                            error('Method not implemented')
                    end
                    
                    % Matching in space
                    dx=xq(1,2)-xq(1,1);
                    dy=yq(2,1)-yq(1,1);
                    ind_x=round((x_tmp/1000-xq(1,1))./dx)+1;
                    ind_y=round((y_tmp/1000-yq(1,1))./dy)+1;
                    
                    ind_xy=(ind_x-1)*size(xq,1)+ind_y;
                    [ind_xy,ind_tmp,ind_ind]=unique(ind_xy);
                    
                    if(compute_deformation)
                        x_obs =[out_xy_tricorner(ind_k(ind_l(ind_m)),1,1);out_xy_tricorner(ind_k(ind_l(ind_m)),2,1);out_xy_tricorner(ind_k(ind_l(ind_m)),3,1)];
                        y_obs =[out_xy_tricorner(ind_k(ind_l(ind_m)),1,2);out_xy_tricorner(ind_k(ind_l(ind_m)),2,2);out_xy_tricorner(ind_k(ind_l(ind_m)),3,2)];
                    else
                        x_obs =[out_xy(ind_k(ind_l(ind_m)),1)];
                        y_obs =[out_xy(ind_k(ind_l(ind_m)),2)];
                    end
                    
                    ind_selected=[];
                    
                    for n=1:size(xq,1),
                        ind_band_rgps=find((y_obs>yq(n,1)-2*res).*(y_obs<=yq(n,1)+2*res));
                        if(length(ind_band_rgps)>0)
                            min_y=min(y_obs(ind_band_rgps));
                            max_y=max(y_obs(ind_band_rgps));
                            min_x=min(x_obs(ind_band_rgps));
                            max_x=max(x_obs(ind_band_rgps));
                            
                            ind_selected=[ind_selected;find((x_tmp/1000>min_x).*(x_tmp/1000<max_x).*(y_tmp/1000>min_y).*(y_tmp/1000<=max_y))];
                        end
                    end
                    
                    ind_selected=unique(ind_selected);
                    nb_selected_point=length(ind_selected);
                    
                    % Compute the deformation in the same way as for the
                    % observations
                    if(compute_deformation && nb_selected_point>min_size_stream)
                        
                        % prepare data for def_from_drift
                        txy2=[x_tmp(ind_selected),y_tmp(ind_selected)]; % m
                        trux2=u_tmp(ind_selected); %m/day
                        trvy2=v_tmp(ind_selected); %m/day
                        
                        % clone of the obs:
                        tdeltat2 =deltats(m);
                        tdnum2   =dnums(l);
                        stream(1:nb_selected_point)=streams(s);
                        
                        % parameters to reject bad triangles after the triangulation
                        min_area=5*1e6;
                        max_area=400*1e6;
                        min_minang=5;
                        max_long_side=25*1e3;
                        showplot=0;
                        filter_noise=0;
                        min_def=0.02;
                        max_level=3;
                        
                        [ntri,ttdnum,ttdeltat,area,xy_tricorner,uv,uv_tricorner,dudx,dudy,dvdx,dvdy,id_stream,quality_index]=def_from_drift(txy2,tdeltat2,tdnum2,trux2,trvy2,streams(s),min_area,max_area,min_minang,max_long_side,showplot,filter_noise,min_def,max_level);
                        
                        if(ntri>0)
                            mobs_data.dnum=[mobs_data.dnum;ttdnum*ones(size(area))];
                            mobs_data.uv = [mobs_data.uv;uv]; % in m/day
                            mobs_data.xy_tricorner=[mobs_data.xy_tricorner;xy_tricorner]; % in km
                            mobs_data.uv_tricorner=[mobs_data.uv_tricorner;uv_tricorner]; % in m/day
                            mobs_data.stream=[mobs_data.stream;id_stream];
                            mobs_data.deltat=[mobs_data.deltat;ttdeltat*ones(size(area))];
                            mobs_data.area=[mobs_data.area;area]; % km^2
                            mobs_data.dudx=[mobs_data.dudx;dudx]; % in /day
                            mobs_data.dudy=[mobs_data.dudy;dudy]; % in /day
                            mobs_data.dvdx=[mobs_data.dvdx;dvdx]; % in /day
                            mobs_data.dvdy=[mobs_data.dvdy;dvdy]; % in /day
                        end
                    else
                        % prepare data for def_from_drift
                        txy2=[x_tmp(ind_selected),y_tmp(ind_selected)]; % m
                        trux2=u_tmp(ind_selected); %m/day
                        trvy2=v_tmp(ind_selected); %m/day
                        
                        % clone of the obs:
                        tdeltat2 =deltats(m);
                        tdnum2   =dnums(l);
                        stream(1:nb_selected_point)=streams(s);
                        
                        
                        mobs_data.x     = [mobs_data.x      ;txy2(:,1)];
                        mobs_data.y     = [mobs_data.y      ;txy2(:,2)];
                        mobs_data.dnum  = [mobs_data.dnum   ;tdnum2*ones(nb_selected_point,1)];
                        mobs_data.uv    = [mobs_data.uv     ;[trux2,trvy2]]; % in m/day
                        mobs_data.stream= [mobs_data.stream ;char(streams(s)*ones(nb_selected_point,1))];
                        mobs_data.deltat= [mobs_data.deltat ;tdeltat2*ones(nb_selected_point,1)];
                    end
                    
                end     % end for each deltat
            end     % end for each date
        end     % end for each stream
    end     % end for dataset
    
    if(compute_deformation)
        if(~isempty(mobs_data.xy_tricorner))
            mobs_data.x = mean(mobs_data.xy_tricorner(:,:,1),2);
            mobs_data.y = mean(mobs_data.xy_tricorner(:,:,2),2);
            mobs_data.scale=sqrt(mobs_data.area);
        end
    end
    
     [ebox_mobs]=regridding(order,xmin,ymin,min_box_size,mobs_data.x/1000,mobs_data.y/1000,mobs_data.stream,mobs_data.dnum-mobs_data.deltat/2,mobs_data.deltat,(mod_date_start+mod_date_end)/2);
end

end