function [moment_1,moment_2, moment_3, tail_exponent]=plots_def_distrib(defo_vec,title_defo)

nb_defo=length(defo_vec);
if(length(title_defo)~=nb_defo)
    error('wrong title_defo')
end

%fig_1=figure;
fig_1bis=figure(1000);
%fig_2=figure;
fig_2bis=figure(1001);
%fig_3=figure;

%ColorOrder=get(gca,'ColorOrder');
ColorOrder=['k';'b';'r';'g'];

for defo_nb=1:nb_defo
    defo=load_defo(defo_vec{defo_nb});
    
    bin_selected=length(defo.data_bin);%-3;
    
    shear=defo.data_bin(bin_selected).invar.shear;
    area=defo.data_bin(bin_selected).area;
    
    length_scale=round(sqrt(mean(area))/1000);
    
    logbins = logspace(-2.5,1,60);
    plotbins = logbins(1:(end-1))+diff(logbins)/2;
    
    dist_obs = histcounts(shear, logbins, 'Normalization', 'pdf');
    
    
    
    figure(fig_2bis)
    subaxis(1,2,1, 'Spacing', 0.06, 'Padding', 0.0, 'Margin', 0.06 );
    loglog(plotbins, dist_obs, '-','Color',[0.2 0.2 0.2]+(i-1)/6.*[0.5,0.5,0.5],'color',ColorOrder(defo_nb))
    hold on
    axis([0.01 0.5 1e-3 1e2])
    set(gca,'fontsize',12);
    
    % Fits
    min_def_fit=0.1;
    max_def_fit=0.5;
    [~,indx0] = min(abs(plotbins-min_def_fit));
    [~,indx1] = min(abs(plotbins-max_def_fit));
    hold on
    x = log10(plotbins(indx0:indx1));
    y = log10(dist_obs(indx0:indx1));
    x(isnan(y))=[];
    y(isnan(y))=[];
    x(isinf(y))=[];
    y(isinf(y))=[];
    C = fit(x',y','poly1')
    xplot = log10([min_def_fit max_def_fit]);
    loglog(10.^xplot, 10.^feval(C, xplot), ['--' ColorOrder(defo_nb)])
    cint = confint(C);
    
    tail_exponent(defo_nb)=    C.p1;
    
    subaxis(1,2,2, 'Spacing', 0.06, 'Padding', 0.0, 'Margin', 0.06 );
    semilogy(plotbins, dist_obs, '-','Color',[0.2 0.2 0.2]+(i-1)/6.*[0.5,0.5,0.5],'color',ColorOrder(defo_nb))
    hold on
    axis([0.01 0.5 1e-3 1e2])
    set(gca,'fontsize',12);
    
    % legend
    xbound=xlim;
    ybound=ylim;
    x_fact=0.5;
    y_fact=0.7/1.5^(defo_nb-1);
    text(x_fact*xbound(2)+(1-x_fact)*xbound(1),y_fact*ybound(2)+(1-y_fact)*ybound(1),[title_defo{defo_nb} ' at ' num2str(length_scale) ' km '],'color',ColorOrder(defo_nb,:),'fontsize',16)
    
    moment_1(defo_nb)=mean(abs(shear));
    moment_2(defo_nb)=mean((abs(shear).^2));
    moment_3(defo_nb)=mean((abs(shear).^3));
    
end