
function [slope_shr,curvature_shr,error_slope_shr]=plots_def_multiscale(defo_vec,title_defo,flags,domain_plot)
% example:
% plots_def_multiscale({'defo_mobs_test_733453.5_733463.5_rgps.mat','defo_733453.5_733463.5_rgps.mat'},{'Mobs','Obs'},[1,1],0,1)

nb_defo=length(defo_vec);
if(length(title_defo)~=nb_defo)
    error('wrong title_defo')
end  

fig_0=figure(1000001);
fig_1=figure(1000002);
fig_2=figure(1000003);
fig_3=figure(1000004);


ColorOrder=get(gca,'ColorOrder');
%ColorOrder=['k';'b';'r';'g']

colors={[1 0.33 1],[1 0 0],[1 0.5 0],[1 1 0],[0.33 1 0],[0 1 1],[0 0 1],[0.33 0 1],[0.5 0 1],[1 0 1]};

nb_bins_selected=6; % for the maps to be ploted
nb_lines=2;



for defo_nb=1:nb_defo
    %defo=defo_vec(defo_nb);
    disp(defo_vec{defo_nb})
    defo=load_defo(defo_vec{defo_nb}); 

    nb_bins=length(defo.data_bin);
    
    bins=nb_bins:-1:nb_bins-nb_bins_selected+1;
    bins2=nb_bins:-1:nb_bins-7; % use the 7 smallest bins only
    nb_columns=(nb_bins_selected+rem(nb_bins_selected,nb_lines))/nb_lines;
    
    if(flags(1))
        %*****************************************************************%
        % Deformation Maps
        %*****************************************************************%    
%         %---------------------------- Total deformation ------------------------------%
%         figure;
%         for bin_nb = bins(1:end);
%             subaxis(nb_lines,nb_columns,nb_bins-bin_nb+1, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.03 ,'MarginRight',0.1);
%             plots_def_map_for_multiscale(defo,'eps',[0 0.08],gray2red,bin_nb,domain_plot);
%             axis_plot=axis;
%             text(mean(axis_plot(1:2)),axis_plot(4)-300,[num2str_round(defo.data_bin_mean.length(bin_nb)/1000) ' km'],'FontSize',16,'Color',colors{bin_nb-1},'FontWeight','bold')
%             rectangle('Position',[axis_plot(1),axis_plot(3),axis_plot(2)-axis_plot(1),axis_plot(4)-axis_plot(3)],...
%                 'LineWidth',4,'EdgeColor',colors{bin_nb-1})
%             axis off
%         end;
%         B=colorbar;
%         set(B, 'Position', [ .9 .33 .02 .37])     
        %---------------------------- Shear ------------------------------%
        figure;
        for bin_nb = bins(1:end);
            subaxis(nb_lines,nb_columns,nb_bins-bin_nb+1, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.03 ,'MarginRight',0.1);
            plots_def_map_for_multiscale(defo,'shear',[0 0.08],gray2red,bin_nb,domain_plot);
            axis_plot=axis;
            text(mean(axis_plot(1:2)),axis_plot(4)-300,[num2str_round(defo.data_bin_mean.length(bin_nb)/1000) ' km'],'FontSize',16,'Color',colors{bin_nb-1},'FontWeight','bold')
            rectangle('Position',[axis_plot(1),axis_plot(3),axis_plot(2)-axis_plot(1),axis_plot(4)-axis_plot(3)],...
                'LineWidth',4,'EdgeColor',colors{bin_nb-1})
            axis off
        end;
        B=colorbar;
        set(B, 'Position', [ .9 .33 .02 .37])     
%         %-------------------------- Divergence ---------------------------%
%         figure;
%         for bin_nb = bins(1:end);
%             subaxis(nb_lines,nb_columns,nb_bins-bin_nb+1, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.03 ,'MarginRight',0.1);
%             plots_def_map_for_multiscale(defo,'div',[-0.04 0.04],blue2red,bin_nb,domain_plot);
%             axis_plot=axis;
%             text(mean(axis_plot(1:2)),axis_plot(4)-300,[num2str_round(defo.data_bin_mean.length(bin_nb)/1000) ' km'],'FontSize',16,'Color',colors{bin_nb-1},'FontWeight','bold')
%             rectangle('Position',[axis_plot(1),axis_plot(3),axis_plot(2)-axis_plot(1),axis_plot(4)-axis_plot(3)],...
%                 'LineWidth',4,'EdgeColor',colors{bin_nb-1})
%             axis off
%         end;
%         B=colorbar;
%         set(B, 'Position', [ .9 .33 .02 .37])    
%         %-------------------------- Vorticity ----------------------------%    
%         figure;
%         for bin_nb = bins(1:end);
%             subaxis(nb_lines,nb_columns,nb_bins-bin_nb+1, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.03 ,'MarginRight',0.1);
%             plots_def_map_for_multiscale(defo,'vor',[-0.06 0.06],blue2red,bin_nb,domain_plot);
%             axis_plot=axis;
%             text(mean(axis_plot(1:2)),axis_plot(4)-300,[num2str_round(defo.data_bin_mean.length(bin_nb)/1000) ' km'],'FontSize',16,'Color',colors{bin_nb-1},'FontWeight','bold')
%             rectangle('Position',[axis_plot(1),axis_plot(3),axis_plot(2)-axis_plot(1),axis_plot(4)-axis_plot(3)],...
%                 'LineWidth',4,'EdgeColor',colors{bin_nb-1})
%             axis off
%         end;
%         B=colorbar;
%         set(B, 'Position', [ .9 .33 .02 .37])
    end;
    
    if(flags(2))
        %*****************************************************************%
        % Spatial scaling
        %*****************************************************************%
%         %------------------------------ Total deformation ----------------------------%
%         [slope_eps_tmp,dif_min_slope_eps,dif_max_slope_eps,polyfit_eps]=plots_def_scaling_for_multiscale(defo,'eps',bins2,multifractal);
        %------------------------------ Shear ----------------------------%
        figure(105);
        [slope_shr_tmp,dif_min_slope_shr,dif_max_slope_shr,H_shr,C1_shr,alpha_shr,polyfit1_shr,polyfit2_shr]=plots_def_scaling_for_multiscale(defo,'shear',bins2,ColorOrder(defo_nb,:));
%         %---------------------------- Divergence -------------------------%
%         [slope_div_tmp,dif_min_slope_div,dif_max_slope_div,polyfit_div]=plots_def_scaling_for_multiscale(defo,'div',bins2,multifractal);
%         %---------------------------- Vorticity ------------------------------%
%         [slope_vor_tmp,dif_min_slope_vor,dif_max_slope_vor,polyfit_vor]=plots_def_scaling_for_multiscale(defo,'vor',bins2,multifractal);

%     slope_eps(defo_nb,:)=slope_eps_tmp;
    mean_shr(defo_nb,:)=slope_shr_tmp;
    slope_shr(defo_nb,:)=slope_shr_tmp;
%     slope_div(defo_nb,:)=slope_div_tmp;
%     slope_vor(defo_nb,:)=slope_vor_tmp;
    
    
%     error_slope_eps(defo_nb,:)=dif_max_slope_eps-dif_min_slope_eps;
    error_slope_shr(defo_nb,:)=dif_max_slope_shr-dif_min_slope_shr;
%     error_slope_div(defo_nb,:)=dif_max_slope_div-dif_min_slope_div;
%     error_slope_vor(defo_nb,:)=dif_max_slope_vor-dif_min_slope_vor;
    
    curvature_shr(defo_nb)=polyfit2_shr(1);
    H_shr
    C1_shr
    alpha_shr

    
        disp('error slope')
%         error_slope_eps(defo_nb,2)
        error_slope_shr(defo_nb,2)
%         error_slope_div(defo_nb,2)
%         error_slope_vor(defo_nb,2)
        
        q=0.5*[0:6];
        x=(q(1):0.01:q(end));
        xbound=[-0.5 3.5];
        ybound=[-0.5 2];
        
%         %------------------------------ Total deformation ----------------------------%
%         figure(fig_0)
%         hold on
%         errorbar(q,-slope_eps_tmp,-dif_max_slope_eps,dif_min_slope_eps,'x','color',ColorOrder(defo_nb,:),'linewidth',1)
%         plot(x,polyfit_eps(1)*x.^2+polyfit_eps(2)*x+polyfit_eps(3),'color',ColorOrder(defo_nb,:),'MarkerSize',16,'linewidth',1)
%         %axis([0 q(end)+0.1 -0.1 3]);
%         ylabel('\beta(q) for eps','fontsize',24);
%         xlabel('q','fontsize',24);
%         set(gca,'fontsize',12);
%         
%         % legend
%         %xbound=xlim;
%         %ybound=ylim;
%         axis([xbound ybound])
%         x_fact=0.05;
%         y_fact=0.95-0.1*(defo_nb-1);
%         text(x_fact*xbound(2)+(1-x_fact)*xbound(1),y_fact*ybound(2)+(1-y_fact)*ybound(1),title_defo{defo_nb},'color',ColorOrder(defo_nb,:),'fontsize',20)
        
        %------------------------------ Shear ----------------------------%
        figure(fig_1)
        hold on
        errorbar(q,-slope_shr_tmp,-dif_max_slope_shr,dif_min_slope_shr,'x','color',ColorOrder(defo_nb,:),'linewidth',1)
        %plot(x,H_shr+,'color',ColorOrder(defo_nb,:),'MarkerSize',16,'linewidth',1)
        %plot(x,polyfit1_shr(1)*x+polyfit1_shr(2),'color',ColorOrder(defo_nb,:),'MarkerSize',16,'linewidth',1)
        plot(x,polyfit2_shr(1)*x.^2+polyfit2_shr(2)*x+polyfit2_shr(3),'color',ColorOrder(defo_nb,:),'MarkerSize',16,'linewidth',1)
        %axis([0 q(end)+0.1 -0.1 3]);
        ylabel('\beta(q) for shr','fontsize',24);
        xlabel('q','fontsize',24);
        set(gca,'fontsize',12);
        
        % legend
        %xbound=xlim;
        %ybound=ylim;
        axis([xbound ybound])
        x_fact=0.05;
        y_fact=0.95-0.1*(defo_nb-1);
        text(x_fact*xbound(2)+(1-x_fact)*xbound(1),y_fact*ybound(2)+(1-y_fact)*ybound(1),title_defo{defo_nb},'color',ColorOrder(defo_nb,:),'fontsize',20)
        
%         %---------------------------- Divergence -------------------------%
%         figure(fig_2)
%         hold on
%         errorbar(q,-slope_div_tmp,-dif_max_slope_div,dif_min_slope_div,'x','color',ColorOrder(defo_nb,:),'linewidth',1)
%         plot(x,polyfit_div(1)*x.^2+polyfit_div(2)*x+polyfit_div(3),'color',ColorOrder(defo_nb,:),'MarkerSize',16,'linewidth',1)
%         %axis([0 q(end)+0.1 -0.1 3]);
%         ylabel('\beta(q) for |div|','fontsize',24);
%         xlabel('q','fontsize',24);
%         set(gca,'fontsize',12);
%         
%         % legend
%         %xbound=xlim;
%         %ybound=ylim;
%         axis([xbound ybound])
%         x_fact=0.05;
%         y_fact=0.95-0.1*(defo_nb-1);
%         text(x_fact*xbound(2)+(1-x_fact)*xbound(1),y_fact*ybound(2)+(1-y_fact)*ybound(1),title_defo{defo_nb},'color',ColorOrder(defo_nb,:),'fontsize',20)
%         
%         %---------------------------- Vorticity ------------------------------%
%         figure(fig_3)
%         hold on
%         errorbar(q,-slope_vor_tmp,-dif_max_slope_vor,dif_min_slope_vor,'x','color',ColorOrder(defo_nb,:),'linewidth',1)
%         plot(x,polyfit_vor(1)*x.^2+polyfit_vor(2)*x+polyfit_vor(3),'color',ColorOrder(defo_nb,:),'MarkerSize',16,'linewidth',1)
%         %axis([0 q(end)+0.1 -0.1 3]);
%         ylabel('\beta(q) for |vor|','fontsize',24);
%         xlabel('q','fontsize',24);
%         set(gca,'fontsize',12);
%         
%         % legend
%         %xbound=xlim;
%         %ybound=ylim;
%         axis([xbound ybound])
%         x_fact=0.05;
%         y_fact=0.95-0.1*(defo_nb-1);
%         text(x_fact*xbound(2)+(1-x_fact)*xbound(1),y_fact*ybound(2)+(1-y_fact)*ybound(1),title_defo{defo_nb},'color',ColorOrder(defo_nb,:),'fontsize',20)
%         
%         %---------------------------- opening/closing ------------------------------%
%         [invar]=invariants(defo.data.dudx,defo.data.dudy,defo.data.dvdx,defo.data.dvdy);
%         size(invar.div)
%         size(defo.data.area)
%         size(defo.data.deltat)
%         divergence_tmp=invar.div.*defo.data.area.*1e-6.*defo.data.deltat;
%         divergence_tmp=divergence_tmp(defo.data.indices);
%         divergence(defo_nb)=sum(divergence_tmp);
%         opening(defo_nb)=sum(divergence_tmp(divergence_tmp()>0));
%         ridging(defo_nb)=sum(divergence_tmp(divergence_tmp<0));
%         total_area(defo_nb)=sum(defo.data.area).*1e-6;
    
    end;
    
    if(flags(3))
        %*****************************************************************%
        % Deformation Maps
        %*****************************************************************%    

        %---------------------------- Shear ------------------------------%
        figure(10+defo_nb);
        bin_nb = bins(1);
        plots_def_map_for_multiscale(defo,'shear',[0.005 0.1],gray2red,bin_nb,domain_plot);
        axis_plot=axis;
        text(mean(axis_plot(1:2)),axis_plot(4)-200,title_defo(defo_nb),'FontSize',16,'Color','k','FontWeight','bold')
        
        axis off
        B=colorbar;
        set(B, 'Position', [ .9 .33 .02 .37])     
%         %-------------------------- Divergence ---------------------------%
%         figure(200+defo_nb);
%         bin_nb = bins(1);
%         plots_def_map_for_multiscale(defo,'div',[-0.05 0.05],blue2red,bin_nb,domain_plot);
%         axis_plot=axis;
%         text(mean(axis_plot(1:2)),axis_plot(4)-200,title_defo(defo_nb),'FontSize',16,'Color','k','FontWeight','bold')
%         axis off
%         B=colorbar;
%         set(B, 'Position', [ .9 .33 .02 .37])    
%         %-------------------------- Vorticity ----------------------------%    
%         figure;
%         for bin_nb = bins(1:end);
%             subaxis(nb_lines,nb_columns,nb_bins-bin_nb+1, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.03 ,'MarginRight',0.1);
%             plots_def_map_for_multiscale(defo,'vor',[-0.06 0.06],blue2red,bin_nb,domain_plot);
%             axis_plot=axis;
%             text(mean(axis_plot(1:2)),axis_plot(4)-300,[num2str_round(defo.data_bin_mean.length(bin_nb)/1000) ' km'],'FontSize',16,'Color',colors{bin_nb-1},'FontWeight','bold')
%             rectangle('Position',[axis_plot(1),axis_plot(3),axis_plot(2)-axis_plot(1),axis_plot(4)-axis_plot(3)],...
%                 'LineWidth',4,'EdgeColor',colors{bin_nb-1})
%             axis off
%         end;
%         B=colorbar;
%         set(B, 'Position', [ .9 .33 .02 .37])
    end;
    
end;

end

%-------------------------------------------------------------------------%
% Subfunctions
%-------------------------------------------------------------------------%

function [str]=num2str_round(scale)
scale=round(scale);
deg=floor(log10(scale));
str=num2str(round((scale/10^(-1+deg)))/10^(1-deg));
end






