%Created by Philipp to compare grids before and after regriding
%Prev and post output needed, can be enabled by setting simul_in.regrid_output_flag to 1
%Old grid in red, new overlayed in black
%Circles show location of lowest angle in old grid.
%If the circles are not visible, zoom in and out once. That normally causes matlab to redraw them correctly
%Expect to zoom in to see anything.
%Example: plot_mesh_regrid('test2',22)
%Example: plot_mesh_regrid('test1',1)

function plot_mesh_regrid(filename,step)

fprev=strcat('regrid_',filename,'_prev',int2str(step))
fprev=sprintf('%s.mat',fprev);
fpost=strcat('regrid_',filename,'_post',int2str(step))
%fpost=strcat('regrid_',filename,'_test')
fpost=sprintf('%s.mat',fpost);


load(fprev);

[mesh_prev] = importbamg(simul_out.bamg.mesh, simul_out.bamg.geom);

node_x_prev = mesh_prev.node.x' + simul_out.UM(1:2:end)*1e-3;
node_y_prev = mesh_prev.node.y' + simul_out.UM(2:2:end)*1e-3;

x_tricorner_prev(:,:) = node_x_prev(mesh_prev.element.num_node);%*1000;
y_tricorner_prev(:,:) = node_y_prev(mesh_prev.element.num_node);%*1000;
indices_prev= 1:length(x_tricorner_prev);


%finding the element that tiggered the regridding happened


boundarynodes = ismember(mesh_prev.element.num_node(:,:),mesh_prev.boundary.from_msh(:,1:2));
boundarysum   = sum(boundarynodes,2);
notlandlocked = find(3>boundarysum);

xy_ang(:,:,1) = node_x_prev(mesh_prev.element.num_node(notlandlocked,:))*1000;
xy_ang(:,:,2) = node_y_prev(mesh_prev.element.num_node(notlandlocked,:))*1000;


[~,minang,~,~,~] = area_min_angle_mex(length(xy_ang),xy_ang);

mn = min(minang);

loc = find(minang == mn)




clear simul_out


if exist(fpost, 'file')
    
    post_flag=1;

    load(fpost);
    
    [mesh_post] = importbamg(simul_out.bamg.mesh, simul_out.bamg.geom);
    
    node_x_post = mesh_post.node.x' ;%+ simul_out.UM(1:2:end)*1e-3;
    node_y_post = mesh_post.node.y' ;%+ simul_out.UM(2:2:end)*1e-3;
    
    x_tricorner_post(:,:) = node_x_post(mesh_post.element.num_node);%*1000;
    y_tricorner_post(:,:) = node_y_post(mesh_post.element.num_node);%*1000;
    indices_post= 1:length(x_tricorner_post);

else
    post_flag=0;

end



figure()
hold on;


h1 =patch(x_tricorner_prev(indices_prev,:)',y_tricorner_prev(indices_prev,:)',zeros(size(x_tricorner_prev))','EdgeColor','red','LineWidth',2,'FaceAlpha',0.0) 

if post_flag==1
  h2 =patch(x_tricorner_post(indices_post,:)',y_tricorner_post(indices_post,:)',zeros(size(x_tricorner_post))','EdgeColor','black','LineWidth',2,'FaceAlpha',0.0) 
end
% plot(mesh_prev.element.x(notlandlocked(loc)),mesh_prev.element.y(notlandlocked(loc)),'ob', 'MarkerSize',25)
% plot(mesh_prev.element.x(notlandlocked(loc)),mesh_prev.element.y(notlandlocked(loc)),'ob', 'MarkerSize',75)
% plot(mesh_prev.element.x(notlandlocked(loc)),mesh_prev.element.y(notlandlocked(loc)),'ow', 'MarkerSize',50)
% plot(mesh_prev.element.x(notlandlocked(loc)),mesh_prev.element.y(notlandlocked(loc)),'ow', 'MarkerSize',5)
% plot(mesh_prev.element.x(notlandlocked(loc)),mesh_prev.element.y(notlandlocked(loc)),'or', 'MarkerSize',100)
% plot(mesh_prev.element.x(notlandlocked(loc)),mesh_prev.element.y(notlandlocked(loc)),'or', 'MarkerSize',10)
% plot(mesh_prev.element.x(notlandlocked(loc)),mesh_prev.element.y(notlandlocked(loc)),'or', 'MarkerSize',1)

axis equal

end

