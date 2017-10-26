function [mesh,element]=msh2mat(filename)
% Read a .msh file and create mesh and element

fid=fopen(filename);

% get the format
fgetl(fid); 
MeshFormat=fgetl(fid);
fgetl(fid);


if(strcmp(MeshFormat,'2.2 0 8'))
    disp(['MeshFormat ', MeshFormat]);
    node_info_nb=4;
    element_info_nb=5;
elseif(strcmp(MeshFormat,'3 0 8'))
    disp(['MeshFormat ', MeshFormat]);
    node_info_nb=5;
    element_info_nb=4;
    for i=1:14
        fgetl(fid);
    end
else
    error(['MeshFormat ', MeshFormat, 'is not implemented'])
end
   
% number of nodes
fgetl(fid);
mesh.Nn=fscanf(fid, ['%d'],1)

% load the nodes
tc=fscanf(fid, ['%f'],node_info_nb*mesh.Nn);
mesh.node.num = tc(1:node_info_nb:end-1);
x   = tc(2:node_info_nb:end-2)';
y   = tc(3:node_info_nb:end-1)';
z   = tc(4:node_info_nb:end)';

% x,y,z -> lat lon
[lat,lon,radius] = cart2lla(x,y,z);
mesh.node.lat=rad2deg(lat);
mesh.node.lon=rad2deg(lon);

%polar stereographic projection
lat=mesh.node.lat;
lon=mesh.node.lon;

m_proj('Stereographic','lon',-45,'lat',90,'radius',60);
[x,y]=m_ll2xy(lon,lat);
mesh.node.x=x*6378.273; 
mesh.node.y=y*6378.273; % in km from the axes centered on the pole

% load the elements and boundary condition (both are objects for Gmsh)
fgetl(fid); fgetl(fid); fgetl(fid);
nb_objects=fscanf(fid, ['%d'],1);

% initialization to the maximum size (nb_objects)
element.num_node=zeros(nb_objects,3);       
mesh.boundary.from_msh=zeros(nb_objects,3);

% Loop over the objects
ind_bc=1;
ind_el=1;
for i=1:nb_objects

    info=fscanf(fid, ['%d'],element_info_nb);
    element_type=info(2)+1;

	if(element_type==1)
        val=fscanf(fid,'%d',element_type);
	elseif(element_type==2)
        val=fscanf(fid,'%d',element_type);
		mesh.boundary.from_msh(ind_bc,1)=val(1);
		mesh.boundary.from_msh(ind_bc,2)=val(2);
        mesh.boundary.from_msh(ind_bc,3)=info(element_info_nb-1);
		ind_bc=ind_bc+1;
	elseif(element_type==3)
        val=fscanf(fid,'%d',element_type);
		element.num_node(ind_el,1)=val(1);
		element.num_node(ind_el,2)=val(2);
		element.num_node(ind_el,3)=val(3);
		ind_el=ind_el+1;
	else
		disp(['reading error: unknown element type? element_type: ',num2str(element_type), ' element id: ',num2str(ind_el)]);
		return;
    end
end

% reduction of the size of element.num_node and mesh.boundary.from_msh
element.num_node      =element.num_node      (1:ind_el-1,:);
mesh.boundary.from_msh=mesh.boundary.from_msh(1:ind_bc-1,:);

% number of elements
mesh.Ne=ind_el-1;

% clean the mesh and recompute the x, y coordinates
[mesh,element]=clean_mesh(mesh,element);

%polar stereographic projection
lat=mesh.node.lat;
lon=mesh.node.lon;

m_proj('Stereographic','lon',-45,'lat',90,'radius',60);
[x,y]=m_ll2xy(lon,lat);
mesh.node.x=x*6378.273; 
mesh.node.y=y*6378.273; % in km from the axes centered on the pole

disp('loading elements finished')

fileout=[filename(1:end-4) '.mat'];
save(fileout,'element','mesh');

return

