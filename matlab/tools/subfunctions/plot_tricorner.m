function plot_tricorner(defo,variable_name,colormap_name,caxis_range,title_figure,domain,box_bound,masked,visible,figure_format,date,mesh_filename)
%Here are a list of various options which can be set
font_size = 16;            %Sets font size of colorbar and date
pic_quality = '-r600';     %resolution for eps, pdf, tiff, and png
pic_resolution_jpeg = 1200; %resolution for everything else, use 1200 for high
plot_grid           =   0; %If not zero the mesh lines are ploted. If zoomed out only the mesh lines will be visible
background_color    = [1 1 1]; %Set background color, typical is either grey [.7 .7 .7] ot white [1 1 1]
plot_domain_grid    =   1; %When 1 the actual domain baoundaries are plotted, closed with black and open 1

if(isempty(defo))
    return;
elseif(isempty(defo.data))
    return;
end

if(masked==1)
    indices=defo.data.indices;
else    
    indices=1:length(defo.data.xy_tricorner);
end

if(strcmp(variable_name,'div') || strcmp(variable_name,'shear') || strcmp(variable_name,'vor') || strcmp(variable_name,'eps'))
    tmp=invariants(defo.data.dudx(indices),defo.data.dudy(indices),defo.data.dvdx(indices),defo.data.dvdy(indices));
elseif(strcmp(variable_name,'u') || strcmp(variable_name,'v') || strcmp(variable_name,'speed') || strcmp(variable_name,'orientation'))
    tmp.u=defo.data.uv(indices,1)/1000;
    tmp.v=defo.data.uv(indices,2)/1000;
    tmp.speed=sqrt(defo.data.uv(indices,1).^2+defo.data.uv(indices,2).^2)/1000;
    tmp.orientation = atan2(defo.data.uv(indices,1),defo.data.uv(indices,2))/pi*180;
else
    tmp_data=getfield(defo.data,variable_name);
    tmp=defo.data;
    tmp=setfield(tmp,variable_name,tmp_data(indices));
end

if(visible)
fig1 =    figure();
else
fig1 =    figure('visible','off');
end

hold on

%Setting background color
whitebg(background_color); %Phil
set(gcf,'InvertHardcopy','off'); %Makes the background be included when the file is saved

tmp_field=getfield(tmp,variable_name)';




%plot(x_coast,y_coast,'k','linewidth',0.5);

%if ~strcmp(domain,'square')
%    m_proj('Stereographic','lon',-45,'lat',90,'radius',60);
%    m_coast('l','patch',[.7 .7 .7]);
%end

if(~isempty(date) && date~=0)
    date_defo=defo.data.dnum(indices)'-0.5*defo.data.deltat(indices)';
    diff_date=date-date_defo;
    diff_date=[diff_date;diff_date;diff_date];
    x_tricorner=(defo.data.xy_tricorner(indices,:,1)'+defo.data.uv_tricorner(indices,:,1)'.*diff_date)/1000;
    y_tricorner=(defo.data.xy_tricorner(indices,:,2)'+defo.data.uv_tricorner(indices,:,2)'.*diff_date)/1000;
else
    x_tricorner=defo.data.xy_tricorner(indices,:,1)'/1000;
    y_tricorner=defo.data.xy_tricorner(indices,:,2)'/1000;
end

if plot_grid == 0
    patchplot = patch(x_tricorner,y_tricorner,tmp_field,'EdgeColor','none');
else 
    patchplot = patch(x_tricorner,y_tricorner,tmp_field,'EdgeColor','green'); %Phil
end
% Coastlines
axis_tmp=axis;
hold on; plot_coastlines_and_boundaries_c(mesh_filename);
axis(axis_tmp)


% % test to highligt selected triangles
% eps=getfield(tmp,'eps')';
% deformed=find(eps>0.02);
% field=getfield(tmp,variable_name)';
% patch(defo.data.xy_tricorner(indices(deformed),:,1)'/1000,defo.data.xy_tricorner(indices(deformed),:,2)'/1000,field(deformed),'EdgeColor','black')

colormap(colormap_name)
if strcmp('pink',colormap_name) & strcmp(domain,'landfast')
    cm=colormap; cm=flipud(cm); colormap(cm);
end


caxis(caxis_range);



box on;

set(gca,'Layer','top')
set(gca,'XTickLabel',{})
set(gca,'YTickLabel',{})
set(gca, 'LooseInset', get(gca, 'TightInset'));

if ischar(domain)
    if strcmp(domain,'kara')
        axis([700 1300 700 1300]);
    elseif strcmp(domain,'bigkara')
        axis([0 3000 -600 1400]);
    elseif strcmp(domain,'bk')
        axis([550 2950 -800 1400]);
    elseif strcmp(domain,'karagate')
        axis([2000 2300 350 650]);
    elseif strcmp(domain,'tipns')
        axis([1050 1550 350 850]);
    elseif strcmp(domain,'beaufort')
        axis([-2450 300 -750 2000]);
    elseif strcmp(domain,'tinykara')
        axis([400,1400,450,1400]);
    elseif strcmp(domain,'minykara')
        axis([935,1010,1175,1250]);
    elseif strcmp(domain,'arctic')
        axis([-2450 1050 -1000 2500]);
    elseif strcmp(domain,'bigarctic')
       axis([-2450 2950 -2900 2500]); 
    elseif strcmp(domain,'square')
       axis([-46 46 -46 46]); 
    elseif strcmp(domain,'squarebig')
       axis([-460 460 -460 460]); 
    elseif strcmp(domain,'squaresmall')
       axis([-46 46 -46 46]); 
    elseif strcmp(domain,'squaresmalltop')
       axis([-15 15 16 46]);  
    elseif strcmp(domain,'MITgcm')
       axis([-3700 3800 -5200 3750]);  
    elseif strcmp(domain,'landfast')
       axis([600 1800 700 1300]);  
    elseif strcmp(domain,'4MIT')
       axis([-1300 350 -2500 -678]);
       %axis equal tight
    elseif strcmp(domain,'BKF')
       axis(1.0e+03*[-0.4    1.3   -2.5    0.3]);       
%       axis equal tight
    else
        axis equal tight
    end
else
    axis(domain)
end

Tick_spacing=500;
axis_tmp=axis;
set(gca,'XTick',[floor(axis_tmp(1)/Tick_spacing):1:ceil(axis_tmp(2)/Tick_spacing)]*Tick_spacing)
set(gca,'YTick',[floor(axis_tmp(3)/Tick_spacing):1:ceil(axis_tmp(4)/Tick_spacing)]*Tick_spacing)

set(gca,'DataAspectRatio',[1 1 1])

 

if(~isempty(box_bound))
    rectangle('Position',[box_bound(1),box_bound(3),box_bound(2)-box_bound(1),box_bound(4)-box_bound(3)],...
     'LineWidth',2,'EdgeColor','g')
end

%Extra boxes
%box_bound = [1050 1550 350 850];
%rectangle('Position',[box_bound(1),box_bound(3),box_bound(2)-box_bound(1),box_bound(4)-box_bound(3)],...
% 'LineWidth',2,'EdgeColor','g');
%box_bound = [700 1300 700 1300];
%rectangle('Position',[box_bound(1),box_bound(3),box_bound(2)-box_bound(1),box_bound(4)-box_bound(3)],...
% 'LineWidth',2,'EdgeColor','g');
%box_bound  = [2000 2300 350 650];
%rectangle('Position',[box_bound(1),box_bound(3),box_bound(2)-box_bound(1),box_bound(4)-box_bound(3)],...
% 'LineWidth',2,'EdgeColor','g');
%
%clear box_bound

allText   = findall(gcf, 'type', 'text');
allAxes   = findall(gcf, 'type', 'axes');
allFont   = [allText; allAxes];
set(allFont,'FontSize',16);

fig1.Units = 'centimeters';
fig1.Position = [0 0 20 20];




axisPos = get(gca, 'OuterPosition');

paperPos = get(gcf, 'PaperPosition');
figureUnits = get(gcf, 'Units');
set(gcf, 'Units', get(gcf,'PaperUnits'));
figurePos = get(gcf, 'Position');




myaxis = axis;
%axis_ratio=axisPos(4)/axisPos(3);
axis_ratio=(myaxis(4)-myaxis(3))/(myaxis(2)-myaxis(1));
figure_ratio=figurePos(4)/figurePos(3);

if(figure_ratio>axis_ratio)
    figurePos_4=figurePos(3)*axis_ratio;
    figurePos_3=figurePos(3);
else
    figurePos_3=figurePos(4)/axis_ratio;
    figurePos_4=figurePos(4);
end

set(gcf, 'PaperSize', [figurePos_3 figurePos_4+0.03]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figurePos_3 figurePos_4]);


set(gcf, 'Position', [figurePos(1) figurePos(2) figurePos_3 figurePos_4]);


if ischar(domain)
    if strcmp(domain,'4MIT') || strcmp(domain,'BKF')
        gcc=colorbar('peer',gca,'horiz','position',[0.6  0.05 0.3   0.05 ]);%,'fontsize',font_size);  
        Xtick_range = [caxis_range(1) 0.5*caxis_range(1)+0.5*caxis_range(2) caxis_range(2)];
        set(gcc,'XTick',Xtick_range,'fontsize',font_size);
        if variable_name == 'c'
            cb_text = 'concentration';
        elseif strcmp(variable_name,'thick')
            cb_text = 'thickness [m]' ;        
        else
            cb_text = variable_name;
        end
        text(0.6+0.05, 0.05+0.075,cb_text,'units','normalized','BackgroundColor',background_color,'FontSize',font_size,'EdgeColor',background_color);
        
    else
        gcc=colorbar('peer',gca,'horiz','position',[0.67 0.91 0.226 0.03 ]);
        set(gcc,'XTick',caxis_range,'fontsize',font_size/1.8)
    end
else
    gcc=colorbar('peer',gca,'horiz','position',[0.67 0.91 0.226 0.03 ]);
    set(gcc,'XTick',caxis_range,'fontsize',font_size/1.8)
end



refresh(gcf)



textstring=datestr(date);
%Adds date
if exist('textstring')==1
    if ischar(domain)
    if strcmp(domain,'4MIT') || strcmp(domain,'BKF')
            text(0.025, 0.1,textstring,'units','normalized','BackgroundColor','white','FontSize',font_size,'EdgeColor','k')
        else
            text(0.05, 0.95,textstring,'units','normalized','BackgroundColor','white','FontSize',font_size,'EdgeColor','k')
        end
    end
else
    text(0.05, 0.95,textstring,'units','normalized','BackgroundColor','white','FontSize',font_size,'EdgeColor','k')
    %annotation(fig1,'textbox',[0.12 0.98 0.01 0.01], 'String',textstring,'BackgroundColor','w', 'EdgeColor','k', 'FitBoxToText','on','VerticalAlignment','top','FontSize',12)
end


%Adding moorings if mooring_data is in defo_temp
%strategy is to check the current_time in the moorings.output to find the
%right one

if isfield(defo,'moorings') 
    
    
    if strcmp(colormap_name,'jet') | strcmp(colormap_name,'parula') 
        arrow_color='k';
    elseif strcmp(colormap_name,'rev_gris')
        arrow_color='r';
    else
        arrow_color='green'
    end

    
    %add a reference arrow of fixed length. 
    rel_x = 0.8;
    rel_y = 0.2;
    ref_arrow.u=100;
    ref_arrow.v=0;
    
    %relative coordinates are calculated
    tempx = get(gca,'xlim');
    tempy = get(gca,'ylim');
    ref_arrow.x=(1-rel_x)*tempx(1)+rel_x*tempx(2);
    ref_arrow.y=(1-rel_y)*tempy(1)+rel_y*tempy(2);
    
    %add to original defo
    defo.moorings.x(end+1) = ref_arrow.x;
    defo.moorings.y(end+1) = ref_arrow.y;
    defo.moorings.VT(end+1)=ref_arrow.u;
    defo.moorings.VT(end+1)=ref_arrow.v;
    
    %Mask to not plot zero values, both x and y values most be zero
    speed     = abs(defo.moorings.VT(1:2:end)) + abs(defo.moorings.VT(2:2:end));
    non_zero_mask = find(speed>0);
    
    
    xlimits = get(gca,'xlim');
    ylimits = get(gca,'ylim');
    
    %text(0.65, 0.2,'20 cm/s        ','units','normalized','BackgroundColor','white','FontSize',font_size,'EdgeColor','k')   
    text(rel_x-0.15, rel_y,'20 cm/s:','units','normalized','BackgroundColor',background_color,'FontSize',font_size,'EdgeColor',background_color)   
    %q = quiver(defo.moorings.x,defo.moorings.y,defo.moorings.VT(1:2:end),defo.moorings.VT(2:2:end),'AutoScale','off','Color',arrow_color)   ;
    q = quiver(defo.moorings.x(non_zero_mask),defo.moorings.y(non_zero_mask),defo.moorings.VT(non_zero_mask*2-1),defo.moorings.VT(non_zero_mask*2),'AutoScale','off','Color',arrow_color)   ;
    
    xlim(xlimits)
    ylim(ylimits)
    
end


    
if(strcmp(figure_format,'eps'))
    print_format='epsc2'
else
    print_format=figure_format;
end


filename=[title_figure((double(title_figure)>=48).*(double(title_figure)<=122)==1),'_', variable_name,'.' figure_format];

if(strcmp(figure_format,'eps')||strcmp(figure_format,'pdf')||strcmp(figure_format,'tiff')||strcmp(figure_format,'png'))
    print(['-d',print_format],pic_quality,'-painters',filename)
    %exportfig(filename,'Resolution',1200,'Renderer','painters','Format',figure_format);
    disp(['Saved file ', filename])
   
elseif(~strcmp(figure_format,''))
    exportfig(filename,'Resolution',pic_resolution_jpeg,'Renderer','painters','Format',figure_format);
    disp(['Saved file ', filename])
end

if(~isempty(box_bound))
   % % to select an area and to plot the histogram of the variables in this specific area
   figure
   bin_length=(caxis_range(2)-caxis_range(1))/40;
   bin_length=0.01;
   hist(tmp_field((defo.data.x>box_bound(1)).*(defo.data.x<box_bound(2)).*(defo.data.y>box_bound(3)).*(defo.data.y<box_bound(4))==1),[caxis_range(1)-bin_length/2:bin_length:caxis_range(2)+bin_length/2])
end

end
