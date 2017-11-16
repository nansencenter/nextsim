
function plot_mesh(filename,domain)
%function plot_mesh(filename,domain, resol)

f=sprintf('%s.mat',filename);
load(f);

if(~exist('lon_coast'))
    load arctic_coasts.mat;
    %load coast_hires_xy.mat;
end

%load([domain resol '.mat'])

node_x=mesh.node.x;
node_y=mesh.node.y;

h=trisurf(element.num_node,node_x,node_y,zeros(size(node_x)),0,'EdgeColor','black'); view(2);
colormap(gray)

hold on;
plot(x_coast,y_coast,'k','linewidth',0.5);
if strcmp(domain,'kara')
    axis([1000 2500 0 1500]);  
elseif strcmp(domain,'bigkara')
    axis([600 3000 -800 1500]);
elseif strcmp(domain,'arctic')
    axis([-2500 1500 -1250 2500]);
else
   axis([-2500 2900 -3000 3000]);
end;
