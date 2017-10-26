
function plots_def_map_for_multiscale(defo,invar_name,caxis_invar,colormap_invar,bin_nb,domain)

if(isempty(defo))
    return;
elseif(isempty(defo.data))
    return;
end

%load arctic_coasts_light.mat;

ls=defo.data_bin(bin_nb).scale/1000;
x_center=defo.data_bin(bin_nb).x/1000;
y_center=defo.data_bin(bin_nb).y/1000;
x_patch=[x_center-ls/1.75,x_center+ls/1.75,x_center+ls/1.75,x_center-ls/1.75]';
y_patch=[y_center-ls/1.75,y_center-ls/1.75,y_center+ls/1.75,y_center+ls/1.75]';
invar=getfield(defo.data_bin(bin_nb).invar,invar_name);
patch(x_patch,y_patch,invar','FaceColor','flat','EdgeColor','none');
hold on;
box on;
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);

plot_coastlines_and_boundaries_c('small_Arctic_10km.msh',false)
%plot(x_coast+5,y_coast-5,'k','linewidth',0.5);
if strcmp(domain,'kara')
    axis([1000 2500 0 1500]);
elseif strcmp(domain,'bigkara')
    axis([600 2900 -800 1500]);
elseif strcmp(domain,'arctic')
    axis([-2500 1500 -1250 2500]);
elseif strcmp(domain,'bigarctic')
    axis([-2500 2900 -2000 3000]);
elseif strcmp(domain,'square')
   axis([-46 46 -46 46]); 
end;
colormap(colormap_invar)
caxis(caxis_invar)
%text(-2100,1400,[num2str_round(defo.data_bin_mean.length(bin_nb)/1000) ' km'],'FontSize',16,'FontWeight','bold')
%title(invar_name)
end

%-------------------------------------------------------------------------%
% Subfunctions
%-------------------------------------------------------------------------%

    function [str]=num2str_round(scale)
        scale=round(scale);
        deg=floor(log10(scale));
        str=num2str(round((scale/10^(-1+deg)))/10^(1-deg));
    end