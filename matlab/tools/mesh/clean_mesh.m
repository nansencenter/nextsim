function [mesh,element] = clean_mesh(mesh,element)
% This function find nodes which are not part of the mesh and remove them
% without "fucking up" the numbering of the nodes.

disp('searching for non connected nodes')

element_num_node_old=([element.num_node(:,1);element.num_node(:,2);element.num_node(:,3)]);
[C,IA,IC]=unique([element.num_node(:,1);element.num_node(:,2);element.num_node(:,3)]);

%Useless nodes are removed
mesh.node.lat=mesh.node.lat(C);
mesh.node.lon=mesh.node.lon(C);

% mesh.node.num is updated
mesh.node.num=C;

%Node numbering of the elements are changed so the numbering
%corresponds to the updated mesh.node
element_num_node_new=IC;
element.num_node(:,1)=element_num_node_new(          1:  mesh.Ne);
element.num_node(:,2)=element_num_node_new(  mesh.Ne+1:2*mesh.Ne);
element.num_node(:,3)=element_num_node_new(2*mesh.Ne+1:3*mesh.Ne);

%Adjusting the node numbering in boundary_msh file to account for the removed nodes.
new_boundary_from_msh=mesh.boundary.from_msh;
for i=1:length(mesh.boundary.from_msh(:,1)),
    for j=1:2,
        new_indice=find(mesh.boundary.from_msh(i,j)==C);
        if(isempty(new_indice))
            new_boundary_from_msh(i,j)=0;
        else
            new_boundary_from_msh(i,j)=new_indice;
        end
    end
end

%Removing unneeded boundary from boundary_msh
unconnected_boundary=find(~(new_boundary_from_msh(:,1).*new_boundary_from_msh(:,2)));
new_boundary_from_msh(unconnected_boundary,:)=[];

mesh.boundary.from_msh=new_boundary_from_msh;

% redefine the node ids for the nodes being on the boundaries to correspond
% to BAMG convention where boundary nodes are the first
old_boundary_from_msh=mesh.boundary.from_msh;
new_boundary_from_msh=mesh.boundary.from_msh;

old_element_num_node=element.num_node;
new_element_num_node=element.num_node;

new_mesh_node_lat=mesh.node.lat;
new_mesh_node_lon=mesh.node.lon;

[old_node_id]=unique([mesh.boundary.from_msh(:,1);mesh.boundary.from_msh(:,2)]);

old_mesh_node_id=(1:length(mesh.node.lat));
new_mesh_node_id=(1:length(mesh.node.lat));

new_mesh_node_lat(old_node_id)=[];
new_mesh_node_lat=[mesh.node.lat(old_node_id),new_mesh_node_lat];
new_mesh_node_lon(old_node_id)=[];
new_mesh_node_lon=[mesh.node.lon(old_node_id),new_mesh_node_lon];
new_mesh_node_id(old_node_id)=[];
new_mesh_node_id=[old_mesh_node_id(old_node_id),new_mesh_node_id];

% % initial solution was much too slow
% for i=1:length(new_mesh_node_id)
%     for j=1:2,
%         new_boundary_from_msh(old_boundary_from_msh(:,j)==new_mesh_node_id(i),j)=old_mesh_node_id(i);
%     end
%     new_element_num_node(old_element_num_node==new_mesh_node_id(i))=old_mesh_node_id(i);
% end

% % % new method, much more rapid and giving the same results
[sorted_new_mesh_node,ind_sort]=sort(new_mesh_node_id);
sorted_old_mesh_node=old_mesh_node_id(ind_sort);
new_boundary_from_msh(:,1:2)=sorted_old_mesh_node(old_boundary_from_msh(:,1:2));
new_element_num_node(:,1:3)=sorted_old_mesh_node(old_element_num_node(:,1:3));

element.num_node=new_element_num_node;
mesh.boundary.from_msh=new_boundary_from_msh;

mesh.node.lat=new_mesh_node_lat;
mesh.node.lon=new_mesh_node_lon;

% Diagnostic
old_Nn  = mesh.Nn;
mesh.Nn = length(C);
disp(['number of non connected nodes removed: ' int2str(old_Nn-mesh.Nn)])

end

