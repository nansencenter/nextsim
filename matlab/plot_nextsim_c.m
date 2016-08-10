function resplot(field,step,domain,dir)

% clearvars -except step;
%field='Velocity';
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

if nargin==3, dir=''; end

%Here are a list of various options which can be set
plot_grid           = 0;            % If not zero the mesh lines are ploted. If zoomed out only the mesh lines will be visible
plot_coastlines     = 1;            % When 1 the actual domain baoundaries are plotted, closed with black and opened in green
plot_date           = 0;            % 0 by default, 1 if we want to display the date on the figure
font_size           = 16;           % Sets font size of colorbar and date
background_color    = [1 1 1];      % white color by default
save_figure         = 1;            % 0 (default) we do not save the figure
figure_format       = 'png';        % can be pdf, tiff, png or jpeg
pic_quality         = '-r600';      % Resolution for eps, pdf, tiff, and png
visible             = 1;            % we display the figure on the screen (may be set to 0 when generating a large amount of figures)

for p=0:0

  [mesh_out,data_out] = neXtSIM_bin_revert(dir, p, step);

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
  
  %we generate the landmask from the grid informations
  if strcmp(domain,'TOPAZ')
    reduced_domain=1;
  else
    reduced_domain=0;
  end;
 
  %To treat the case of plotting the velocity (2D) field
  i=1; %always valid when the plotting a 1D field
  if strcmp(field,'M_VTu')
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

  %--------------------------------------------------------------------------------------------------------
  % We arrange the figure in an "optimal" manner using subfunctions (appended at the bottom of this script)
  %--------------------------------------------------------------------------------------------------------
  %
  set_region_adjustment(domain);
  %
  set_axis_colormap_colorbar(domain,field,v);
  %
  set_figure_cosmetics(font_size,background_color);
  
  %We can now color the ocean in blue...
  patch(x(:,mask_water)/1000,y(:,mask_water)/1000,[0 0.021 0.53],'EdgeColor','none');
  
  %We plot the coastlines and boundaries
  if plot_coastlines == 1
      plot_coastlines_and_boundaries(domain);
  end;
  
  %------------------------------
  % We save the figure (optional)
  %------------------------------
  if save_figure
      filename = sprintf('neXtSIM_%s_%d',field,step);
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
    if (strcmp(domain,'TOPAZ'))
    load topazreducedsplit2.mat
    elseif (strcmp(domain,'MITgcm_4km'))
    load MITgcm4kmsplit2.mat
    elseif (strcmp(domain,'MITgcm_9km'))
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
    
    [closed_boundaryX,closed_boundaryY]= mapll(closed_boundaryLat,closed_boundaryLon,60,-45,'N');
    [free_boundaryX,free_boundaryY]= mapll(free_boundaryLat,free_boundaryLon,60,-45,'N');
    %Plotting closed mesh boundaries (e.g. coastlines)
    plot(closed_boundaryX,closed_boundaryY,'Color',[0.3 0.3 0.3],'LineWidth',0.2);
    %Plotting open mesh boundaries
    plot(free_boundaryX,free_boundaryY,'c','LineWidth',2);
end
function set_axis_colormap_colorbar(domain,field,v)
    %We set the axis limits and the colormap
    if (strcmp(field,'Concentration'))
        caxis([0 1]);
        load('ice_conc_cmap64.mat')
        colormap(ice_conc_cmap64);
    elseif (strcmp(field,'Thickness'))
        caxis([0, 4]);
        colormap('parula');
    elseif strcmp(field,'Damage')
        caxis([0.995, 1]);
        load('ice_damage_cmap128.mat')
        colormap(ice_damage_cmap128);
    elseif strcmp(field,'M_VT')
        caxis([-0.10, 0.10]);
        colormap('blue2red');
    elseif strcmp(field,'Divergence')
        caxis([min(min(v{:})), max(max(v{:}))]);
        colormap('blue2red');
    elseif strcmp(field,'Snow')
        caxis([0, 0.5]);
        colormap('parula');
    end;

    % We predefine the position of the colorbar (only made for TOPAZ and MITgcm) that fits well the
    % figure and we add it to the figure
    if strcmp(domain,'TOPAZ')
        colorbar_Xposition=0.715;
        colorbar_Yposition=0.65;
        width_scale_factor=0.5;
        height_scale_factor=0.3;
        % We add a colorbar (settings below are (at least) good for MITgcm or Topaz setups)
        c=colorbar;
        cpos = c.Position;
        cpos(1) = colorbar_Xposition;
        cpos(2) = colorbar_Yposition;
        cpos(3) = width_scale_factor*cpos(3);
        cpos(4) = height_scale_factor*cpos(4);
        c.Position=cpos;
    elseif strcmp(domain,'MITgcm')
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
    else
        % We add a colorbar (settings below are (at least) good for MITgcm or Topaz setups)
        c=colorbar;
    end;
end    
function set_region_adjustment(domain)
        if strcmp(domain,'kara')
            axis([700 1300 700 1300]);
        elseif strcmp(domain,'arch')
            %axis equal tight
            axis([0 170 -20 290])
            %axis([0 160 -20 290])
        elseif strcmp(domain,'bigkara')
            axis([600 2900 -800 1500]);
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
            axis([-3570 3720 -6100 3750]);
        elseif strcmp(domain,'landfast')
            axis([600 1800 700 1300]);
        elseif strcmp(domain,'4MIT')
            axis([-1300 350 -2500 -678]);
        elseif strcmp(domain,'BKF')
            axis equal tight
        elseif strcmp(domain,'FramStrait')
            axis([200 1200 -1400 0]);
        elseif strcmp(domain,'TOPAZ')
            axis([-2800 3800 -4800 2250]);
        end;
end
function set_figure_cosmetics(font_size,background_color)

        whitebg(background_color);
        set(gcf,'InvertHardcopy','off'); %Makes the background be included when the file is saved

        box on;
        set(gca,'Layer','top')
        set(gca,'XTickLabel',{})
        set(gca,'YTickLabel',{})
        %set(gca, 'LooseInset', get(gca, 'TightInset'));

        Tick_spacing=300;
        axis_tmp=axis;
        set(gca,'XTick',[floor(axis_tmp(1)/Tick_spacing):1:ceil(axis_tmp(2)/Tick_spacing)]*Tick_spacing)
        set(gca,'YTick',[floor(axis_tmp(3)/Tick_spacing):1:ceil(axis_tmp(4)/Tick_spacing)]*Tick_spacing)

        set(gca,'DataAspectRatio',[1 1 1],'LineWidth',0.3)

        allText   = findall(gcf, 'type', 'text');
        allAxes   = findall(gcf, 'type', 'axes');
        allFont   = [allText; allAxes];
        set(allFont,'FontSize',font_size);
end

%%%% The following does not work, but we keep it for possible later use %%%

%[xmask,ymask,pmask]=load_landmask(domain,reduced_domain);
  %pmask(pmask==0)=NaN;
function [xmask,ymask,pmask]=load_landmask(domain,reduced_domain)

if strcmp(domain,'TOPAZ')
        % open the topography file of TOPAZ
        try
            fid=fopen('regional.grid.b');
            A=fgetl(fid); idm=sscanf(A,'%f');
            A=fgetl(fid); jdm=sscanf(A,'%f');
            fclose(fid);
        catch
            error('regional.grid.b cannot be found');
            return
        end
        try
            plon=loada('regional.grid.a',1,idm,jdm);
            plat=loada('regional.grid.a',2,idm,jdm);
            pdepth=loada('regional.depth.a',1,idm,jdm);
            pmask=double(pdepth~=pdepth(1,1));
        catch
            error('regional.depth.a and/or regional.grid.a not found');
        end
        
        [xmask,ymask]=mapll(plat,plon,90,-45,'N');
        
        % Buldozer operations
        %We remove the few remaining grid cells lying the Mediterranean Sea
        f=find(plat<=50 & plon>=0);
        pmask(f)=0;
        
        % Bering Strait is closed in TOPAZ4. To have the opportunity to open it,
        % the mask of the cells beyond Bering Strait are set to 2
        pmask(222,813:836)=2;
        
        % we select the indices i  and j of the subdomain if
        % domain_reduced=1
        if(reduced_domain)
            %We select the domain that will be meshed
            min_i=200;
            max_i=min(630,idm);
            min_j=435;         
            max_j=jdm;
            pmask=pmask(min_i:max_i,min_j:max_j);
            xmask=xmask(min_i:max_i,min_j:max_j);
            ymask=ymask(min_i:max_i,min_j:max_j);
        end;
        
elseif strcmp(domain,'MITgcm_9km')
        try
            pmask =readbin('9km_hFacC.data',[420*2 384*2]);
        catch
            disp('9km_hFacC.data cannot be not found');
            return
        end
        
        if(reduced_domain)
            % We select the domain that will be meshed
            min_i=240;
            max_i=min(470,size(pmask,1));
            min_j=1;
            max_j=min(200,size(pmask,2));
            pmask=pmask(min_i:max_i,min_j:max_j);
        end
               
elseif strcmp(domain,'MITgcm_4km')
        try
            pmask =readbin('hFacC.data',[420*4 384*4]);
        catch
            disp('hFacC.data cannot be found');
            return
        end
        
        if(reduced_domain)
            % We select the domain that will be meshed
            % Phil: Reduced domain here correspond to the neXtSIMF forecast domain
            min_i=480;
            min_i=320;
            max_i=min(940,size(pmask,1));
            min_j=1;
            max_j=min(400,size(pmask,2));
            pmask=pmask(min_i:max_i,min_j:max_j);
        end
end;

end


