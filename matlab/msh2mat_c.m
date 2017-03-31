function [mesh,element]=msh2mat_c(filename)
% Read a .msh file and create mesh and element

fid=fopen(filename);
fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid);
fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid);

% number of nodes
mesh.Nn=fscanf(fid, '%d',1);

% load the nodes
tc=fscanf(fid, '%f',4*mesh.Nn);
mesh.node.num = tc(1:4:end-1);
x   = tc(2:4:end-2)';
y   = tc(3:4:end-1)';
z   = tc(4:4:end)';

% x,y,z -> lat lon
[lat,lon] = cart2lla(x,y,z);
mesh.node.lat=rad2deg(lat);
mesh.node.lon=rad2deg(lon);

% load the elements and boundary condition (both are objects for Gmsh)
fgetl(fid); fgetl(fid); fgetl(fid);
nb_elements=fscanf(fid, '%d',1);

% initialization to the maximum size (nb_objects)
element.num_node=zeros(nb_elements,3);       
mesh.boundary.from_msh=zeros(nb_elements,3);

% Loop over the objects
ind_bc=1;
ind_el=1;
for i=1:nb_elements

    info=fscanf(fid, '%d',5);
    element_type=info(2)+1;

    if(element_type==1)
        val=fscanf(fid,'%d',element_type);
    elseif(element_type==2)
        val=fscanf(fid,'%d',element_type);
        mesh.boundary.from_msh(ind_bc,1)=val(1);
        mesh.boundary.from_msh(ind_bc,2)=val(2);
        mesh.boundary.from_msh(ind_bc,3)=info(4);
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

return

