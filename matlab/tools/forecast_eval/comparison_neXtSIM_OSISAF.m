forecast_date_str='2015-10-26';
forecast_name='BKF3km_forecast';

% step to analyse
forecast_step_array=[2:5]; % [2,7] By default from 2 to 7

close all

limit_concentration=0.6; % Not very useful when the quality flag is set to 30

% initialization
mean_error=zeros(1,length(forecast_step_array));
median_error=mean_error;
min_error=mean_error;
max_error=mean_error;

mean_error_persistence=mean_error;
median_error_persistence=mean_error;
min_error_persistence=mean_error;
max_error_persistence=mean_error;

percentage_good_forecast=mean_error;
percentage_good_persistence=mean_error;

% data time interval in days
data_dt=2; % in days

% define path
define_default_data_path;

forecast_dnum=datenum(forecast_date_str);
forecast_date=datestr(forecast_dnum,'yyyymmdd');

% %----- load the forecast simul_in ----- (not needed at the moment)
% filename=['simul_in_', forecast_name, '_', forecast_date, '.mat'];
% full_path=[forecast_dir, forecast_date, '/',filename];
% file=dir(full_path);
% if(isempty(file))
%     disp(full_path)
%     disp('No input files found!')
%     return
% end
% load(full_path);

%-------------------------------------------------------------------------%
% Loading obs data for persistence
%-------------------------------------------------------------------------%
% From documentation:
% "on day 0 around 0600 UTC, low-resolution ice drift datasets are distributed 
% which cover the period from day -3 to day -1"
% When our forecast is available we have then access to the observation for
% forecast_dnum-2 to forecast_dnum

disp('Persistence')
startdate=datestr(forecast_dnum-2,'yyyymmdd')
enddate  =datestr(forecast_dnum,'yyyymmdd')

startyear = startdate(1:4);
startmonth = startdate(5:6);
startday = startdate(7:8);
endyear = enddate(1:4);
endmonth = enddate(5:6);
endday = enddate(7:8);

filename=['ice_drift_nh_polstere-625_multi-oi_', startdate, '1200-', enddate, '1200.nc'];
full_path=[indir 'OSISAF_ice_drift/', endyear,'/', endmonth ,'/',filename];
file=dir(full_path);
if(isempty(file))
    disp(full_path)
    disp('No input files found!')
    return
end

% load displacement and status flag
dX_persistence=ncread(full_path,'dX'); % km
dY_persistence=ncread(full_path,'dY'); % km
status_flag_persistence=ncread(full_path,'status_flag'); % - (30= nominal quality)

% Loop over the steps
for array_i=1:length(forecast_step_array),
    
    figure

    forecast_step=forecast_step_array(array_i);
    
    %----- define information on the date -----
    if(forecast_step<data_dt) % need to be equal or higher than data_dt
        error(['forecast_step (', num2str(forecast_step), ') needs to be equal or higher than data_dt (', num2str(data_dt),')'])
    end
    
    startdate=datestr(forecast_dnum+forecast_step-data_dt,'yyyymmdd')
    enddate  =datestr(forecast_dnum+forecast_step        ,'yyyymmdd')
    
    startyear = startdate(1:4);
    startmonth = startdate(5:6);
    startday = startdate(7:8);
    endyear = enddate(1:4);
    endmonth = enddate(5:6);
    endday = enddate(7:8);
    
    %-------------------------------------------------------------------------%
    % Loading model data
    %-------------------------------------------------------------------------%
    node_id=[];
    node_x=[];
    node_y=[];
    UM=[];
    
    % we load only the first and last
    step_ids=[forecast_step-data_dt,forecast_step];
    for i=1:2
        
        % load the step i
        filename=['simul_out_', forecast_name, '_', forecast_date, '_step', num2str(step_ids(i)), '.mat'];
        full_path=[forecast_dir, forecast_date, '/',filename];
        file=dir(full_path);
        if(isempty(file))
            disp(full_path)
            disp('No input files found!')
            return
        end
        load(full_path);
        
        % Filter as a function of the minimum concentration in the
        % surrounding elements
        node_reduced=(simul_out.ind_reduced(1:2:end-1)+1)/2;
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
    end
    
    % as we loaded only the first and last step:
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
    
    %-------------------------------------------------------------------------%
    % Loading obs data
    %-------------------------------------------------------------------------%
    filename=['ice_drift_nh_polstere-625_multi-oi_', startdate, '1200-', enddate, '1200.nc'];
    full_path=[indir 'OSISAF_ice_drift/', endyear,'/', endmonth ,'/',filename];
    file=dir(full_path);
    if(isempty(file))
        disp(full_path)
        disp('No input files found!')
        return
    end
    
    % load xc and yc
    xc=ncread(full_path,'xc'); % km
    yc=ncread(full_path,'yc'); % km
    
    dx=62.5;
    
    % load displacement and status flag
    dX=ncread(full_path,'dX'); % km
    dY=ncread(full_path,'dY'); % km
    status_flag=ncread(full_path,'status_flag'); % - (30= nominal quality)
    
    % % load lat/lon
    % lat=ncread(filename,'lat'); % deg
    % lon=ncread(filename,'lon'); % deg
    %
    % % Check if the following projection is consistent with xc and yc
    % m_proj('Stereographic','lon',-45,'lat',90,'radius',70);
    % [x,y]=m_ll2xy(lon,lat);
    % x=(x*6378.273)/1.0307; % km
    % y=(y*6378.273)/1.0304; % km
    %
    % tolerance=10; % km
    % max_error_x=max(abs(x(1:end,1)-xc ));
    % max_error_y=max(abs(y(1,1:end)-yc'));
    % max_error=max(max_error_x,max_error_y);
    % if(max_error>tolerance)
    %     error(['The projection used generates an error of ',num2str(max_error),' km'])
    % end
    
    %-------------------------------------------------------------------------%
    % Comparison of the data
    %-------------------------------------------------------------------------%
    
    % position of the model nodes in km and with a correction due to the sligthly different projection
    x_start=(x_start/1000)/1.0307;
    y_start=(y_start/1000)/1.0304;
    
    x_end=(x_end/1000)/1.0307;
    y_end=(y_end/1000)/1.0304;
    
    % Compute the indices on the 62.5 km regular grid of OSISAF
    indice_x=round((x_start-xc(1))/dx)+1;
    indice_y=round((yc(1)-y_start)/dx)+1;
    
     % Compute the indices on the original 12.5 km regular grid of OSISAF
    xc_1_12_5=xc(1)-dx;
    yc_1_12_5=yc(1)+dx;
    dx_12_5=12.5;
     
    indice_x_12_5=round((x_start  -xc_1_12_5)/dx_12_5)+1;
    indice_y_12_5=round((yc_1_12_5-y_start  )/dx_12_5)+1;
    
    % not very efficient to make it this way, but ok for the moment
    error_tmp=[]; 
    error_tmp_persistence=[]; 
    
    % Loop over the regular grid
    for i=1:length(xc),
        for j=1:length(yc),
            % only look at where there is data (we could aslo use the quality flags in addition)
            if(~isnan(dX(i,j)) && ~isnan(dY(i,j)) && (status_flag(i,j)==30))
                %indices=find((indice_x==i).*(indice_y==j));
                
                local_indice_x_12_5=indice_x_12_5-(i-1)*5;
                local_indice_y_12_5=indice_y_12_5-(j-1)*5;
               
                % find the local indices that are in a 11 by 11 pixels area
                % minus the 12 corner pixels (see Lavergne et al. 2010)
                indices=find(...
                    (local_indice_x_12_5>=1 ).*...
                    (local_indice_x_12_5<=11).*...
                    (local_indice_y_12_5>=1 ).*...
                    (local_indice_y_12_5<=11).*...
                    (2~=((local_indice_x_12_5==1 )+(local_indice_y_12_5==1 ))).*...
                    (2~=((local_indice_x_12_5==1 )+(local_indice_y_12_5==2 ))).*...
                    (2~=((local_indice_x_12_5==1 )+(local_indice_y_12_5==10))).*...
                    (2~=((local_indice_x_12_5==1 )+(local_indice_y_12_5==11))).*...
                    (2~=((local_indice_x_12_5==11)+(local_indice_y_12_5==1 ))).*...
                    (2~=((local_indice_x_12_5==11)+(local_indice_y_12_5==2 ))).*...
                    (2~=((local_indice_x_12_5==11)+(local_indice_y_12_5==10 ))).*...
                    (2~=((local_indice_x_12_5==11)+(local_indice_y_12_5==11))).*...
                    (2~=((local_indice_x_12_5==2 )+(local_indice_y_12_5==1 ))).*...
                    (2~=((local_indice_x_12_5==2 )+(local_indice_y_12_5==11))).*...
                    (2~=((local_indice_x_12_5==10)+(local_indice_y_12_5==1 ))).*...
                    (2~=((local_indice_x_12_5==10)+(local_indice_y_12_5==11))));
                  
                % only look at where we have more than 25 model nodes
                if(length(indices)>25)
                    mean_x_start=mean(x_start(indices));
                    mean_y_start=mean(y_start(indices));
                    mean_x_end=mean(x_end(indices));
                    mean_y_end=mean(y_end(indices));
                    dX_mean=mean_x_end-mean_x_start;
                    dY_mean=mean_y_end-mean_y_start;
                    
                    dist_mean_start_c=hypot(mean_x_start-xc(i),mean_x_start-xc(i));
                    
                    % only look at where the mean initial position are
                    % as close as 6.25 km
                    if(dist_mean_start_c<dx/10);
                        hold on
                        
%                         % plot each individual node motion (just for checking)
%                         plot(x_start(indices),y_start(indices),'.')
%                         dX_nodes=x_end(indices)-x_start(indices);
%                         dY_nodes=y_end(indices)-y_start(indices);
%                         %quiver(x_start(indices),y_start(indices),dX_nodes,dY_nodes,0,'Color','k')
                        
                        quiver(mean_x_start,mean_y_start,dX_mean,dY_mean,0,'Color','b')
                                                
                        quiver(xc(i),yc(j),dX(i,j),dY(i,j),0,'Color','r')
                        error_tmp=[error_tmp,hypot(dX_mean-dX(i,j),dY_mean-dY(i,j))];
                        
                        if(~isnan(dX_persistence(i,j)) && ~isnan(dY_persistence(i,j)) && (status_flag_persistence(i,j)==30))
                            quiver(xc(i),yc(j),dX_persistence(i,j),dY_persistence(i,j),0,'Color','g');
                            error_tmp_persistence=[error_tmp_persistence,hypot(dX_persistence(i,j)-dX(i,j),dY_persistence(i,j)-dY(i,j))];
                        end
%                         pause
                    end
                end
            end
        end
    end
    
    mean_error(1,array_i)=mean(error_tmp);
    median_error(1,array_i)=median(error_tmp);
    max_error(1,array_i)=max(error_tmp);
    min_error(1,array_i)=min(error_tmp);
    
    mean_error_persistence(1,array_i)=mean(error_tmp_persistence);
    median_error_persistence(1,array_i)=median(error_tmp_persistence);
    max_error_persistence(1,array_i)=max(error_tmp_persistence);
    min_error_persistence(1,array_i)=min(error_tmp_persistence);
    
    accuracy_OSISAF = 2*4; % 2 standard deviations
    percentage_good_forecast(1,array_i)=sum(error_tmp<accuracy_OSISAF)/length(error_tmp);
    percentage_good_persistence(1,array_i)=sum(error_tmp_persistence<accuracy_OSISAF)/length(error_tmp_persistence);
    
    hold on
    if(~exist('x_coast'))
        load arctic_coasts_light.mat;
    end
    plot(x_coast,y_coast,'k','linewidth',0.5);
    
    axis([0 1500 0 1500]);
    title(['Simulated (blue), observed (red) and persistent 48-hours displacement for step ', num2str(forecast_step-data_dt), ' to ',num2str(forecast_step)])
    
end

% plot evolution of the mean/max/min error
figure
hold on
title('Median/max error')
%plot(forecast_step_array,mean_error,'k:')
plot(forecast_step_array,median_error,'b')
plot(forecast_step_array,max_error,'b--')
plot(forecast_step_array,median_error_persistence,'g')
plot(forecast_step_array,max_error_persistence,'g--')
legend('Median error','Max error','Persistance median error','Persistance max error')
axis([0,7,0,65])
xlabel('step (days)')
ylabel('Displacement error (km)')

% plot evolution of the mean/max/min error
figure
hold on
title(['Success rate (error smaller than ', num2str(accuracy_OSISAF) ,' km (2 OSISAF std))'])
%plot(forecast_step_array,mean_error,'k:')
plot(forecast_step_array,percentage_good_forecast*100,'b')
plot(forecast_step_array,percentage_good_persistence*100,'g')
legend('Forecast success rate','Persistance  success rate')
axis([0,7,0,100])
xlabel('step (days)')
ylabel('Success rate (%)')
