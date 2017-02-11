function plot_nextsim_c(field,step,domain,region_of_zoom,is_sequential,dir,save_figure)

% example of usage:
% plot_nextsim_c('Concentration',4,'topaz',[],true)
% plot_nextsim_c('Damage',4,'topaz',[],true)
% plot_nextsim_c('M_VT',4,'topaz',[],true)

% clearvars -except step;
%field='M_VT';
%field='mld';
%field='Concentration';
%field='Thickness';
%field='SSS';
%field='SST';
%field='Wind';
%field='Ocean';
%field='Vair_factor';
%field='Voce_factor';
%field='Damage';
%field='bathy';

if nargin==5, dir=''; end

%Here are a list of various options which can be set
plot_grid           = 0;            % If not zero the mesh lines are ploted. If zoomed out only the mesh lines will be visible
plot_coastlines     = 1;            % When 1 the actual domain baoundaries are plotted, closed in light gray and opened in cyan. 
                                    % Note though that plotting the coastlines or the grid makes the figure much heavier
plot_date           = 1;            % 0 by default, 1 if we want to display the date on the figure
font_size           = 14;           % Sets font size of colorbar and date
background_color    = [0.85 0.85 0.85];% white color [1 1 1] by default. A substitute could be gray [0.5 0.5 0.5]
figure_format       = 'png';        % can be pdf, tiff, png or jpeg
pic_quality         = '-r300';      % Resolution for eps, pdf, tiff, and png
visible             = 1;            % we display the figure on the screen (may be set to 0 when generating a large amount of figures)

if nargin<7, save_figure = 0; end;     % 0 (default) we do not save the figure

for p=0:0

  if(is_sequential)
      [mesh_out,data_out] = neXtSIM_bin_revert(dir, [], step);
  else
      [mesh_out,data_out] = neXtSIM_bin_revert(dir, p, step);
  end
  %reshape
  var_mx=mesh_out.Nodes_x(mesh_out.Elements);
  var_my=mesh_out.Nodes_y(mesh_out.Elements);
  nr = size(var_mx,1);
  Ne=nr/3;
  Nn=length(mesh_out.Nodes_x);
  x=reshape(var_mx,[3,Ne]);
  y=reshape(var_my,[3,Ne]);
  
  %ice mask and water mask extraction
  mask=data_out.Concentration;
  mask_ice=find(mask>0);
  mask_water=find(mask==0);
  
  i=1;
  %In case we want to plot the velocity or one of its components
  if strcmp(field,'M_VT')
      i=3;
  elseif strcmp(field,'M_VTu')
      field='M_VT';
      i=1;
  elseif strcmp(field,'M_VTv')
      field='M_VT';
      i=2;
  end;
  
  %---------------------------
  % We extract the data fields
  %---------------------------
  field_tmp=data_out.(field);
  length(field_tmp)
  if(length(field_tmp)==Ne)
    v{1}=[field_tmp,field_tmp,field_tmp]';
  elseif(length(field_tmp)==2*Nn)
    var_mc=field_tmp(mesh_out.Elements);
    v{1}=reshape(var_mc,[3,Ne]);
    var_mc=field_tmp(mesh_out.Elements+Nn);
    v{2}=reshape(var_mc,[3,Ne]);
    v{3}=hypot(v{1},v{2});
  elseif(length(field_tmp)==Nn)
    var_mc=field_tmp(mesh_out.Elements);
    v{1}=reshape(var_mc,[3,Ne]);
  else
    error('Not the right dimensions')
  end
  
  %-------------------------------------------------------------------------------------------------------------------
  % We set on a new figure, which we display or not (useful when creating and saving on disk a large bunch of figures)
  %-------------------------------------------------------------------------------------------------------------------
  if(visible)
      fig=figure;
  else
      fig=figure('visible','off');
  end
  
  %-----------------------------------------------------------------------
  % We plot the actual data on the current figure, as well as the landmask
  %-----------------------------------------------------------------------
  if plot_grid == 0
      patch(x(:,mask_ice)/1000,y(:,mask_ice)/1000,v{i}(:,mask_ice),'FaceColor','flat','EdgeColor','none')
  else
      patch(x(:,mask_ice)/1000,y(:,mask_ice)/1000,v{i}(:,mask_ice),'FaceColor','flat','EdgeColor',[0.3 0.3 0.3],'LineWidth',0.2)
  end;
  hold on;

  %----------------------------------------------------------------------------------------------------------------------
  % We arrange the figure in an "optimal" manner using subfunctions (you can check them out at the bottom of this script)
  %----------------------------------------------------------------------------------------------------------------------
  %
  set_region_adjustment(domain,region_of_zoom);
  %
  set_axis_colormap_colorbar(domain,field,v,i,region_of_zoom);
  %
  set_figure_cosmetics(data_out,domain,region_of_zoom,plot_date,background_color,font_size);
  
  %We can now color the ocean in blue...
  patch(x(:,mask_water)/1000,y(:,mask_water)/1000,[0 0.021 0.53],'EdgeColor','none');
  
  %We plot the coastlines and boundaries (optional).
  if plot_coastlines == 1
      plot_coastlines_and_boundaries(domain);
  end;
  
  %------------------------------
  % We save the figure (optional)
  %------------------------------
  if save_figure
      set(fig,'Color',[1 1 1]);
      filename=sprintf('neXtSIM_%s_%d',field,step);
      %Call export_fig to save figure
      if strcmp(figure_format,'png') || strcmp(figure_format,'jpg')
          if isempty(pic_quality)
              pic_quality = '-r300'; %if figure saved as a pdf, this will have no impact
          end;
          export_fig(filename,pic_quality);
      elseif strcmp(figure_format,'pdf')
          export_fig(filename)
      end;
  end;
end;
end

function plot_coastlines_and_boundaries(domain)
% The files defining the domain are on /Data/sim/data/mesh and have to be
% in your matlab path
    if (strcmp(domain,'topaz'))
    load topazreducedsplit2.mat
    elseif (strcmp(domain,'topaz_matthias'))
    load topaz_matthias_split2.mat   
    elseif (strcmp(domain,'mitgcm4km'))
    load MITgcm4kmsplit2.mat
    elseif (strcmp(domain,'mitgcm9km'))
    load MITgcm9kmsplit2.mat
    end
    flag_boundary_fix=10000; %may need to be changed depending on the meshfile (here works for topaz and mit ones)
    flag_boundary_free=10001; %may need to be changed depending on the meshfile (here works for topaz and mit ones)
    boundary   = mesh.boundary.from_msh;
    node_lat   = mesh.node.lat;
    node_lon   = mesh.node.lon;
    %Selecting closed boundaries
    fix = find(flag_boundary_fix==boundary(:,3));
    %Selecting free boundaries
    free   = [];
    for loop_i=1:length(flag_boundary_free)
        fbf = flag_boundary_free(loop_i);
        free= [free;find(fbf==boundary(:,3))];
    end
    closed_boundaryLat  = node_lat(boundary(fix ,1:2,1))';
    closed_boundaryLon  = node_lon(boundary(fix ,1:2,1))';
    free_boundaryLat = node_lat(boundary(free ,1:2,1))';
    free_boundaryLon = node_lon(boundary(free ,1:2,1))';
    
    mppfile = which('NpsNextsim.mpp');

    % projection
    [closed_boundaryX,closed_boundaryY]= mapx_forward(mppfile,closed_boundaryLon(:)',closed_boundaryLat(:)');
    [free_boundaryX,free_boundaryY]= mapx_forward(mppfile,free_boundaryLon(:)',free_boundaryLat(:)');
    
    %reshape and scaling to km
    closed_boundaryX=reshape(0.001*closed_boundaryX,size(closed_boundaryLon));
    closed_boundaryY=reshape(0.001*closed_boundaryY,size(closed_boundaryLon));
    free_boundaryX=reshape(0.001*free_boundaryX,size(free_boundaryLon));
    free_boundaryY=reshape(0.001*free_boundaryY,size(free_boundaryLon));
    
    %Plotting closed mesh boundaries (e.g. coastlines)
    plot(closed_boundaryX,closed_boundaryY,'Color',[0.3 0.3 0.3],'LineWidth',0.2);
    %Plotting open mesh boundaries
    plot(free_boundaryX,free_boundaryY,'g','LineWidth',2);
end
function set_axis_colormap_colorbar(domain,field,v,i,region_of_zoom)
    
    %We set the axis limits, the colormap and set the name for the colorbar
    if (strcmp(field,'Concentration'))
        caxis([0 1]);
        load('ice_conc_cmap64.mat')
        colormap(ice_conc_cmap64);
        name_colorbar='Concentration';
    elseif (strcmp(field,'Thickness'))
        caxis([0, 4]);
        colormap('parula');
        name_colorbar='Thickness (m)';
    elseif strcmp(field,'Damage')
        caxis([0.9995, 1]);
        load('ice_damage_cmap128.mat')
        colormap(ice_damage_cmap128);
        name_colorbar='Damage';
    elseif (strcmp(field,'M_VT') && i==3)
        caxis([0, 0.8]);
        load('ice_speed_cmap128.mat')
        colormap(ice_speed_cmap128);
        name_colorbar='Speed (m/s)';
    elseif (strcmp(field,'M_VT') && i==1)
        caxis([-0.25, 0.25]);
        colormap('blue2red');
        name_colorbar='Speed Ux (m/s)';
    elseif (strcmp(field,'M_VT') && i==2)
        caxis([-0.25, 0.25]);
        colormap('blue2red');
        name_colorbar='Speed Uy (m/s)';
    elseif strcmp(field,'Divergence')
        caxis([min(min(v{:})), max(max(v{:}))]);
        colormap('blue2red');
        name_colorbar='Divergence rate (/day)';
    elseif strcmp(field,'Snow')
        caxis([0, 0.5]);
        colormap('parula');
        name_colorbar='Snow thickness (m)';
    elseif strcmp(field,'Cohesion')
        colormap('parula');
        name_colorbar='Cohesion (Pa)';
    else
        colormap('parula');
        name_colorbar='';
    end;

    % We predefine the position of the colorbar (only made for TOPAZ and MITgcm) that fits well the
    % figure and we add it to the figure
    if (strcmp(domain,'topaz') && isempty(region_of_zoom))
        colorbar_Xposition=0.715;
        colorbar_Yposition=0.63;
        width_scale_factor=0.5;
        height_scale_factor=0.3;
        % We add a colorbar (settings below are (at least) good for MITgcm or Topaz setups)
        c=colorbar;
        cpos = c.Position;
        cpos(1) = colorbar_Xposition;
        cpos(2) = colorbar_Yposition;
        cpos(3) = width_scale_factor*cpos(3);
        cpos(4) = height_scale_factor*cpos(4);
        c.Position = cpos;
        set(get(c,'title'),'string',name_colorbar);
        
    elseif ((strcmp(domain,'mitgcm4km')||strcmp(domain,'mitgcm9km')) && isempty(region_of_zoom))
        colorbar_Xposition=0.67;%(position in the left/right direction. to be modified by hand to your convenience)
        colorbar_Yposition=0.65;
        width_scale_factor=0.5;%(width of the colorbar. to be modified to your convenience)
        height_scale_factor=0.3;%(height of the colorbar. to be modified to your convenience)
        % We add a colorbar
        c=colorbar;
        cpos = c.Position;
        cpos(1) = colorbar_Xposition;
        cpos(2) = colorbar_Yposition;
        cpos(3) = width_scale_factor*cpos(3);
        cpos(4) = height_scale_factor*cpos(4);
        c.Position=cpos;    
        set(get(c,'title'),'string',name_colorbar);
    
    else
        % We add a colorbar (settings below are (at least) good for MITgcm or Topaz setups)
        c=colorbar;
        set(get(c,'title'),'string',name_colorbar);
    end;
end    
function set_region_adjustment(domain,region_of_zoom)
        %%first we adjust depending on the domain
        if strcmp(domain,'topaz')
            axis([-2800 3800 -4800 2250]);
        elseif strcmp(domain,'mitgcm4km')||strcmp(domain,'mitgcm9km')
            axis([-3570 3720 -6100 3750]);
        elseif strcmp(domain,'4MIT')
            axis([-1300 350 -2500 -678]);
        elseif strcmp(domain,'BKF')
            axis equal tight
        elseif strcmp(domain,'kara')
            axis([700 1300 700 1300]);
        elseif strcmp(domain,'bigkara')
            axis([600 2900 -800 1500]);
        elseif strcmp(domain,'bk')
            axis([550 2950 -800 1400]);
        elseif strcmp(domain,'arch')
            %axis equal tight
            axis([0 170 -20 290]);
% Other domains we have been using but that I commented here below for now. Feel free to use though.        
%         elseif strcmp(domain,'tinykara')
%             axis([400,1400,450,1400]);
%         elseif strcmp(domain,'minykara')
%             axis([935,1010,1175,1250]);        
%         elseif strcmp(domain,'square')
%             axis([-46 46 -46 46]);
%         elseif strcmp(domain,'squarebig')
%             axis([-460 460 -460 460]);
%         elseif strcmp(domain,'squaresmall')
%             axis([-46 46 -46 46]);
%         elseif strcmp(domain,'squaresmalltop')
%             axis([-15 15 16 46]); 
        end;
        if ~isempty(region_of_zoom)
            %%Then we adjust depending on chosen region to zoom in
            if strcmp(region_of_zoom,'framstrait')
                axis([200 1200 -1400 0]);
            elseif strcmp(region_of_zoom,'naresstrait')
                axis([-400 100 -950 -450]);
            elseif strcmp(region_of_zoom,'karagate')
                axis([1800 2300 350 650]);
            elseif strcmp(region_of_zoom,'beaufort')
                axis([-2450 -500 -750 1100]);
            elseif strcmp(region_of_zoom,'central_arctic')
                axis([-2450 1050 -1250 2050]);
            elseif strcmp(region_of_zoom,'extended_arctic')
                axis([-2450 2950 -2900 2500]);
            elseif strcmp(region_of_zoom,'kara_landfast')
                axis([600 1800 600 1200]);
            elseif strcmp(region_of_zoom,'tip_novazem')
                axis([1050 1550 250 700]);
            else
                error('This region of zoom does not exist')
            end;
        end;
end
function set_figure_cosmetics(data_out,domain,region_of_zoom,plot_date,background_color,font_size)

        %Adds date
        if plot_date == 1 && isfield(data_out,'Time')
            date=data_out.Time;
            textstring = datestr(date,'yyyy/mm/dd HH:MM');
            if strcmp(domain,'4MIT') || strcmp(domain,'BKF')
                text(0.025, 0.1,textstring,'units','normalized','BackgroundColor','white','FontSize',font_size,'EdgeColor','k')
            elseif strcmp(domain,'arch') 
                text(0.1, 0.03,textstring,'units','normalized','BackgroundColor','white','FontSize',font_size,'EdgeColor','k')
            elseif strcmp(domain,'topaz_matthias') || strcmp(region_of_zoom,'central_arctic')
                text(0.027, 0.95,textstring,'units','normalized','BackgroundColor','white','FontSize',font_size,'EdgeColor','k')
            else
                text(0.2, 0.95,textstring,'units','normalized','BackgroundColor','white','FontSize',font_size,'EdgeColor','k')
            end;
        end;
        whitebg(background_color);
        %set(gcf,'InvertHardcopy','off'); %Makes the background be included when the file is saved
        set(gca,'XColor','k','YColor','k');
        box on;
        %set(gca,'Layer','top')
        set(gca,'XTickLabel',{})
        set(gca,'YTickLabel',{})
        %set(gca, 'LooseInset', get(gca, 'TightInset')); note: was used by philipp for his forecast plots
        Tick_spacing=200;
        axis_tmp=axis;
        set(gca,'XTick',(floor(axis_tmp(1)/Tick_spacing):1:ceil(axis_tmp(2)/Tick_spacing))*Tick_spacing)
        set(gca,'YTick',(floor(axis_tmp(3)/Tick_spacing):1:ceil(axis_tmp(4)/Tick_spacing))*Tick_spacing)
        set(gca,'DataAspectRatio',[1 1 1],'LineWidth',0.3)
        
        %We set the size of the fonts used in all the text present on the
        %figure
        allText   = findall(gcf, 'type', 'text');
        allAxes   = findall(gcf, 'type', 'axes');
        allFont   = [allText; allAxes];
        set(allFont,'FontSize',font_size);
end



