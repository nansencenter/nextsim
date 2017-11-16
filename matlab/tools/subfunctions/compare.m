function [u1,...
    u2,...
    v1,...
    v2,...
    x1,...
    x2,...
    y1,...
    y2,...
    dnum1,...
    dnum2,...
    speed_1,...
    speed_2,...
    i_box,...
    j_box,...
    xmin,...
    ymin,...
    min_box_size,...
    ratio_dom_resol,...
    rmse_u,correlation_u,nb_data_u,...
    rmse_v,correlation_v,nb_data_v]=compare(defo_vec,title_defo,plot_figures)
% example:
% compare({'defo_mobs_test12_rgps.mat','defo_rgps.mat'},{'Model','obs'});

if(length(defo_vec)~=2)
    error('The length of the input vector should be 2')
end

% parameters that could be changed
max_speed_to_plot=40;           % limit of the axis
plot_the_consistency_plots=0;   % set to 1 if you want to check that the temporal and spatial mapping are right

% load of the data
[defo_1,ebox_1]=load_defo(defo_vec{1}); 
[defo_2,ebox_2]=load_defo(defo_vec{2});     

data_1=defo_1.data;
data_2=defo_2.data;

xmin=ebox_1.xmin;
ymin=ebox_1.ymin;
min_box_size=ebox_1.min_box_size;
ratio_dom_resol=ebox_1.ratio_dom_resol;

% selection of the boxes (here we take all the domain)
tmp_ebox_1=ebox_1.mask;
tmp_ebox_2=ebox_2.mask;

% search the corresponding indices
[i_tmp,j_tmp]=find(tmp_ebox_1.*tmp_ebox_2);

% Loop over the data
indices_1=[];
indices_2=[];
i_box=[];
j_box=[];

for i=1:length(i_tmp)

    tmp_indices_1=ebox_1.full{i_tmp(i),j_tmp(i)}';
    tmp_indices_2=ebox_2.full{i_tmp(i),j_tmp(i)}';
    
    out_1_stream=data_1.stream(tmp_indices_1);
    out_1_dnum=data_1.dnum(tmp_indices_1);
    out_1_deltat=data_1.deltat(tmp_indices_1);
    
    out_2_stream=data_2.stream(tmp_indices_2);
    out_2_dnum=data_2.dnum(tmp_indices_2);
    out_2_deltat=data_2.deltat(tmp_indices_2);
    
    % for each stream
    streams=unique(out_1_stream);
    for s=1:length(streams)
        ind_k=find(out_1_stream==streams(s));
        ind_2_k=find(out_2_stream==streams(s));
        
        if(isempty(ind_2_k)) break; end
        
        % for each date
        dnums=unique(out_1_dnum(ind_k));
        for l=1:length(dnums)
            ind_l=find(out_1_dnum(ind_k)==dnums(l));
            ind_2_l=find(out_2_dnum(ind_2_k)==dnums(l));
            
            if(isempty(ind_2_l)) break; end
            
            % for each deltat
            deltats=unique(out_1_deltat(ind_k(ind_l)));
            for m=1:length(deltats)
                ind_m=find(out_1_deltat(ind_k(ind_l))==deltats(m));
                ind_2_m=find(out_2_deltat(ind_2_k(ind_2_l))==deltats(m));
                
                if(isempty(ind_2_m)) break; end
                
                nb_match=min(length(ind_m),length(ind_2_m));
                indices_1=[indices_1;tmp_indices_1(ind_k(ind_l(ind_m(1:nb_match))))];
                indices_2=[indices_2;tmp_indices_2(ind_2_k(ind_2_l(ind_2_m(1:nb_match))))];
                i_box=[i_box,i_tmp(i)*ones(1,nb_match)];
                j_box=[j_box,j_tmp(i)*ones(1,nb_match)];
            end
        end
    end
end

u1=data_1.uv(indices_1,1)/1000; % in km/day
u2=data_2.uv(indices_2,1)/1000; % in km/day
v1=data_1.uv(indices_1,2)/1000; % in km/day
v2=data_2.uv(indices_2,2)/1000; % in km/day
x1=data_1.x(indices_1)/1000; % in km
x2=data_2.x(indices_2)/1000; % in km
y1=data_1.y(indices_1)/1000; % in km
y2=data_2.y(indices_2)/1000; % in km
dnum1=data_2.dnum(indices_2,1);
dnum2=data_2.dnum(indices_2,1);
speed_1=sqrt(u1.^2+v1.^2);
speed_2=sqrt(u2.^2+v2.^2);
error_velocity=sqrt((u2-u1).^2+(v2-v1).^2);
relative_error_velocity=error_velocity./speed_2;

if(plot_figures)

    if(plot_the_consistency_plots)
        figure;
        plot(y1,y2,'.')
        
        figure;
        plot(x1,x2,'.')
        
        figure;
        plot(x1,y1,'.')
        
        figure;
        plot(dnum1,dnum2,'.')
        
        figure;
        plot(double(data_1.stream(indices_1,1)),double(data_2.stream(indices_2,1)),'.')
    end
    
    load arctic_coasts_light.mat;
    
%     name=title_defo{1};
%     plot_data(u1,x1,y1,name,' ice u-velocity (km/day) from ',-max_speed_to_plot,max_speed_to_plot)
%     hold on
%     plot(x_coast,y_coast,'k','linewidth',0.5);
%     
%     plot_data(v1,x1,y1,name,' ice v-velocity (km/day) from ',-max_speed_to_plot,max_speed_to_plot)
%     hold on
%     plot(x_coast,y_coast,'k','linewidth',0.5);
%     
%     plot_data(speed_1,x1,y1,name,' ice speed (km/day) from ',-max_speed_to_plot,max_speed_to_plot)
%     hold on
%     plot(x_coast,y_coast,'k','linewidth',0.5);
%     
%     name=title_defo{2};
%     plot_data(u2,x2,y2,name,' ice u-velocity (km/day) from ',-max_speed_to_plot,max_speed_to_plot)
%     hold on
%     plot(x_coast,y_coast,'k','linewidth',0.5);
%     
%     plot_data(v2,x2,y2,name,' ice v-velocity (km/day) from ',-max_speed_to_plot,max_speed_to_plot)
%     hold on
%     plot(x_coast,y_coast,'k','linewidth',0.5);
%     
%     plot_data(speed_2,x2,y2,name,' ice speed (km/day) from ',-max_speed_to_plot,max_speed_to_plot)
%     hold on
%     plot(x_coast,y_coast,'k','linewidth',0.5);
    
step_vector=100;

diffu=u1-u2;
diffv=v1-v2;

x_ref_vector=-2000;
y_ref_vector=-1500;

s=subplot(1,3,1);
q1=quiver(x2(1:step_vector:end),y2(1:step_vector:end),u1(1:step_vector:end)*20,v1(1:step_vector:end)*20,0,'color','k','AutoScaleFactor',1);
hold on;
set(q1,'MaxHeadSize',0.05);
q2=quiver(x_ref_vector,y_ref_vector,20*10,0,0,'color','k','LineWidth',2);
set(q2,'MaxHeadSize',1e5);
text(x_ref_vector+270,y_ref_vector,'10km/day','FontSize',14)
title('Model','FontSize',18);
%plot(0,0,'ko','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','k')

s=subplot(1,3,2);
q1=quiver(x2(1:step_vector:end),y2(1:step_vector:end),u2(1:step_vector:end)*20,v2(1:step_vector:end)*20,0,'color','k','AutoScaleFactor',1);
hold on;
set(q1,'MaxHeadSize',0.05);
q2=quiver(x_ref_vector,y_ref_vector,20*10,0,0,'color','k','LineWidth',2);
set(q2,'MaxHeadSize',1e5);
text(x_ref_vector+270,y_ref_vector,'10km/day','FontSize',14)
title('Observations','FontSize',18);
%plot(0,0,'ko','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','k')


s=subplot(1,3,3);
q1=quiver(x2(1:step_vector:end),y2(1:step_vector:end),-diffu(1:step_vector:end)*20,-diffv(1:step_vector:end)*20,0,'color','k','AutoScaleFactor',1);
hold on;
set(q1,'MaxHeadSize',0.05);
q2=quiver(x_ref_vector,y_ref_vector,20*10,0,0,'color','k','LineWidth',2);
set(q2,'MaxHeadSize',1e5);
text(x_ref_vector+270,y_ref_vector,'10km/day','FontSize',14)
title('Observations minus Model','FontSize',18);
%plot(0,0,'ko','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','k')


% cumulative probability of the error
figure;
cdfplot(error_velocity)
axis([0 max_speed_to_plot 0 1])

end

% scatter plot for the u component of the velocity
[rmse_u,correlation_u,nb_data_u]=plot_scatter_syl(u2,u1,'Obs','Model',' ice v-velocity (km/day) from ',-max_speed_to_plot,max_speed_to_plot)

% scatter plot for the v component of the velocity
[rmse_v,correlation_v,nb_data_v]=plot_scatter_syl(v2,v1,'Obs','Model',' ice v-velocity (km/day) from ',-max_speed_to_plot,max_speed_to_plot)