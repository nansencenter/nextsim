function [slope,dif_min_slope,dif_max_slope,H,C1,alpha,p1,p2]=plots_def_scaling_for_multiscale(defo,invar_name,bins2,color_tmp)
% plot multifractal spatial scaling

if(isempty(defo))
    return;
elseif(isempty(defo.data))
    return;
end

% bounds for the analysis
%first_bin_for_slope=bins2(5);
%last_bin_for_slope=bins2(2);
first_bin_for_slope=bins2(7);
last_bin_for_slope=bins2(4);


if(strcmp(invar_name,'div'))
    title_name='divergence';
    y_label='Divergence rate';
    small_name='|div|';
elseif(strcmp(invar_name,'shear'))
    title_name='shear';
    y_label='Shear rate';
    small_name='shr';
elseif(strcmp(invar_name,'leads'))
    title_name='lead';
    y_label='Leads ratio';
    small_name='lead';
elseif(strcmp(invar_name,'vor'))
    title_name='vorticity';
    y_label='Vorticity';
    small_name='|vor|';
elseif(strcmp(invar_name,'eps'))
    title_name='total_deformation';
    y_label='Total_deformation';
    small_name='eps';
else
    title_name=invar_name;
    y_label=invar_name;
end

invar_mean = getfield(defo.data_bin_mean,invar_name);

q=0.5*[0:6];
slope=zeros(length(q),1);
min_slope=zeros(length(q),1);
max_slope=zeros(length(q),1);
for j=1:length(q),
    for bin_nb=bins2(1:end)
        if(~isnan(defo.data_bin_mean.length(bin_nb)))
            if(round(q(j))==q(j))
                loglog(defo.data_bin_mean.length(bin_nb)/1000,abs(invar_mean(bin_nb,j)),'o','MarkerEdgeColor','k',...
                    'MarkerFaceColor',color_tmp,'MarkerSize',8,'linewidth',1)
                
                hold on;
            end
        end
    end
    x=log(defo.data_bin_mean.length(first_bin_for_slope:last_bin_for_slope)/1000);
    y=log(abs(invar_mean(first_bin_for_slope:last_bin_for_slope,j)))';
    p=polyfit(x,y,1);
    max_slope(j)=max((y(2:end)-y(1:end-1))./(x(2:end)-x(1:end-1)));
    min_slope(j)=min((y(2:end)-y(1:end-1))./(x(2:end)-x(1:end-1)));
    slope(j)=p(1);
    x=defo.data_bin_mean.length(first_bin_for_slope:last_bin_for_slope)/1000;
    
    if(round(q(j))==q(j))
        loglog(x,exp(p(2))*x.^p(1),'--','linewidth',2,'Color',color_tmp);
    end
    if(q(j)==1)
        disp(['Slope q=' num2str(q(j)) ': ' num2str(p(1))]);
    end
end;
axis([5,1000,1e-6,5e-2]);
ylabel(['<' small_name '^{q}> [day^{-q}]  '],'fontsize',24);
%    ylabel(['' title_name ' rate (day^{-q})  '],'fontsize',24);
xlabel('Spatial scale [km]','fontsize',24);
%set(gca,'XTickLabel',[10 100 1000],'fontsize',12);
box on;

dif_min_slope=min_slope-slope;
dif_max_slope=max_slope-slope;

p1=polyfit(q,-slope',1);
p2=polyfit(q,-slope',2);

curvature=p2(1);
slope_q=p2(2);
disp(['Curvature q: ' num2str(curvature)]);
disp(['slope q: ' num2str(slope_q)]);

box on;

H=-slope(q==1);
C1=(slope(q==1.5)-slope(q==0.5))/(1.5-0.5) + H;
alpha=(slope(q==1.5)-2*slope(q==1)+slope(q==0.5))/(0.5^2)/C1;

end