
function Myconvert_arctic_mesh(filename)
% create mesh and element structures
% If filename ends with .msh, the mesh is loaded from filename
% If filename ends by .a, the mesh is generated from regional.grid.a and regional.depth.a
% If filename ends by .data, the mesh is generated from XC.data,... and hFacC.data
% The new mesh and element are saved in a .mat file, 
% having the same name as filename
% examples: Myconvert_arctic_mesh('MITgcm.data'),
% Myconvert_arctic_mesh('TOPAZ_2_4.a'), Myconvert_arctic_mesh('bigarctic10km.msh')

showplot=1;

disp('Creating a mesh from:')
disp(filename)

if(strcmp(filename(end-3:end),'.msh'))
    % Loading a mesh from gmsh 
    [mesh,element]=msh2mat(filename);
    mesh_name=filename(1:end-4);
elseif(strcmp(filename(end-1:end),'.a')||strcmp(filename(end-1:end),'.b'))
    % Generating the mesh from a regular grid 
    refinement=1;    % by default (refinement=1) we refine the mesh in selected areas (Kara gate and CAA). Only works with split_factor=2
    reduced_domain=1;% by default (reduced_domain=1) we reduce the mesh to the arctic region.
    
    split_factor=2;
    [mesh,element]=grid2mat(filename,showplot,'TOPAZ',refinement,reduced_domain,split_factor);
    mesh_name=filename(1:end-2);
    if(reduced_domain)
        mesh_name=[mesh_name,'reduced'];
    end
    if(refinement)
        mesh_name=[mesh_name,'refined'];     
    else
        mesh_name=[mesh_name,'split' num2str(split_factor)];
    end
elseif(strcmp(filename(end-4:end),'.data'))
    % Generating the mesh from a MITgcm grid 
    refinement=0; % by default we do not refine the mesh but it works.
    reduced_domain=0;% by default (reduced_domain=0) we mesh the whole MITgcm Arctic cube sphere.
    
    % select the MITgcm grid
    if(strcmp(filename(end-7:end-5),'9km'))
        MITgcm_grid='MITgcm_9km'
        split_factor=2;
    else
        MITgcm_grid='MITgcm_4km'
        split_factor=2;
        reduced_domain=1;% by default (reduced_domain=0) we mesh the whole MITgcm Arctic cube sphere.
        %split_factor=4;
        %split_factor=8;
    end
      
    % build the mesh
    [mesh,element]=grid2mat(filename,showplot,MITgcm_grid,refinement,reduced_domain,split_factor);
    mesh_name=filename(1:end-5);
    if(reduced_domain)
        mesh_name=[mesh_name,'reduced'];
    end
    if(refinement)
        mesh_name=[mesh_name,'refined'];     
    else
        mesh_name=[mesh_name,'split' num2str(split_factor)];
    end
elseif(strcmp(filename,'simplesquare'))
    % Generating the mesh for a simple square 
    refinement=0; % by default we do not refine the mesh but it works.
    reduced_domain=0; % by default (reduced_domain=0) we mesh the whole simplesquare domain.
    split_factor=2;
    [mesh,element]=grid2mat(filename,showplot,'simplesquare',refinement,reduced_domain,split_factor);
    mesh_name=filename;
    if(reduced_domain)
        mesh_name=[mesh_name,'reduced'];
    end
    if(refinement)
        mesh_name=[mesh_name,'refined'];     
    else
        mesh_name=[mesh_name,'split' num2str(split_factor)];
    end
elseif(strcmp(filename,'archbox'))
    % Generating the mesh for a simple square 
    refinement=0; % by default we do not refine the mesh but it works.
    reduced_domain=0; % by default (reduced_domain=0) we mesh the whole archbox domain.
    split_factor=2;
    [mesh,element]=grid2mat(filename,showplot,'archbox',refinement,reduced_domain,split_factor);
    mesh_name=filename;
    if(reduced_domain)
        mesh_name=[mesh_name,'reduced'];
    end
    if(refinement)
        mesh_name=[mesh_name,'refined'];     
    else
        mesh_name=[mesh_name,'split' num2str(split_factor)];
    end
else
    error([filename 'is not a valid filename'])
end

% clean the mesh: remove non-connected nodes and boundaries
[mesh,element]=clean_mesh(mesh,element);

%polar stereographic projection
lat=mesh.node.lat;
lon=mesh.node.lon;

m_proj('Stereographic','lon',-45,'lat',90,'radius',60);
[x,y]=m_ll2xy(lon,lat);
mesh.node.x=x*6378.273; 
mesh.node.y=y*6378.273; % in km from the axes centered on the pole

if(showplot)
    node_x=mesh.node.x;
    node_y=mesh.node.y;
    figure
    h=trisurf(element.num_node,node_x,node_y,zeros(size(node_x)),0,'EdgeColor','black'); view(2);hold on
    xbound=xlim;
    ybound=ylim;
    
    ColorOrder=get(gca,'ColorOrder');
    boundary_flags=unique(mesh.boundary.from_msh(:,3));
    for i_flag=1:length(boundary_flags)
        ind_fix=find(mesh.boundary.from_msh(:,3)==boundary_flags(i_flag));
        boundary.num_node=[mesh.boundary.from_msh(ind_fix,1),mesh.boundary.from_msh(ind_fix,2),mesh.boundary.from_msh(ind_fix,1)];
        h=trisurf(boundary.num_node,node_x,node_y,zeros(size(node_x)),0,'EdgeColor',ColorOrder(i_flag,:),'LineWidth',4); view(2); hold on    
        
        x_fact=0.8;
        y_fact=0.95-0.05*(i_flag-1);
        text(x_fact*xbound(2)+(1-x_fact)*xbound(1),y_fact*ybound(2)+(1-y_fact)*ybound(1),['flag: ',num2str(boundary_flags(i_flag))],'color',ColorOrder(i_flag,:),'fontsize',15)
    end
end

if(isnan(mesh.node.x(1)))
	%HACK for meshes with 2 first points badly defined
	mesh.node.x=mesh.node.x(3:end); mesh.node.y=mesh.node.y(3:end);
	mesh.Nn=mesh.Nn-2;
	element.num_node=element.num_node-2; 
	mesh.boundary.from_msh=mesh.boundary.from_msh-2;
end

disp('Compute the shape coef...')
element=shape_coef(mesh,element); % we generate the form parameters of each element    

% Calculate lon and lat
[element.lon,element.lat]=m_xy2ll(element.x/6378.273,element.y/6378.273);

% 
% disp('Building the smoothing matrix... It takes time...')
% % create a smoothing matrix with a radius of 50km
% smoothing_matrix=sparse(mesh.Nn,mesh.Nn);
% l2=50^2;
% x=mesh.node.x;
% y=mesh.node.y;
% for i=1:mesh.Nn,
%     dist2=(x(i)-x).^2+(y(i)-y).^2;
%     ind=(dist2)<l2;
%     weight=(l2-dist2(ind))/l2;
%     smoothing_matrix(i,ind)=weight/sum(weight);
% end

%mesh.smoothing_matrix=smoothing_matrix;
% figure 
% spy(smoothing_mesh)

fileout=sprintf('%s.mat',mesh_name);
save(fileout,'element','mesh');

return
