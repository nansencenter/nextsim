function plots_vel_distrib(defo_vec,title_defo,bins,masked)

nb_defo=length(defo_vec);
if(length(title_defo)~=nb_defo)
    error('wrong title_defo')
end

figure
hold on

ColorOrder=get(gca,'ColorOrder');


for defo_nb=1:nb_defo
    defo=load_defo(defo_vec{defo_nb}); 
    
    if(masked==1)
        indices=defo.data.indices;
    else
        indices=1:length(defo.data.xy_tricorner);
    end
    
    % velocity speed distribution
    speed=0.001*sqrt(defo.data.uv(indices,1).^2+defo.data.uv(indices,2).^2);
    [h,x]=hist(speed,bins);
    plot(x,h/sum(h),'o-','color',ColorOrder(defo_nb,:));
           
    xlim([bins(1) bins(end)])
    
    ylabel('Density function','fontsize',16);
    xlabel('Velocity speed (km/day)','fontsize',16);
    % legend
    xbound=xlim;
    ybound=ylim;
    x_fact=0.5;
    y_fact=0.7-0.1*(defo_nb-1);
    text(x_fact*xbound(2)+(1-x_fact)*xbound(1),y_fact*ybound(2)+(1-y_fact)*ybound(1),title_defo{defo_nb},'color',ColorOrder(defo_nb,:),'fontsize',20)
    
    mean_speed=mean(speed);
    disp(['mean speed rgps: ' num2str(mean_speed)])
end

%to be tried: nhist(data,'binfactor',0.1,'samebins','separate','pdf');