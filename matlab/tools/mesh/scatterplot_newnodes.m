function plot_newnodes(simul_outfile)
%Created by Philipp 2014-11
%Draws a line of the mesh boundary, and a scatter plot of all new nodes
%formed during the last remeshing before the outputfile was generated

load(simul_outfile);

prevnumber = simul_out.bamg.mesh.PreviousNumbering;

%Grabing bamg mesh
[mesh] = importbamg(simul_out.bamg.mesh, simul_out.bamg.geom);


% Adding displacement
node_x = mesh.node.x' + simul_out.UM(1:2:end)*1e-3;
node_y = mesh.node.y' + simul_out.UM(2:2:end)*1e-3;

node_x = node_x./6378.273;
node_y = node_y./6378.273;



%Finding out which nodes are new

%Boundary nodes come first and need to be cut
boundarynodes = length(mesh.boundary.from_msh)+1;

newnodes = find(0==prevnumber);
noborder = find(boundarynodes<newnodes);

newnode_x = node_x(newnodes(noborder));
newnode_y = node_y(newnodes(noborder));



%Plotting
figure();

hold on


%Plotting mesh boundary
boundary = mesh.boundary.from_msh;

%Selecting only closed boundaries
b1 = find(1==boundary(:,3));

%Selecting only open boundaries
b2 = find(0==boundary(:,3));
X = node_x(boundary(b1,1:2,1));
Y = node_y(boundary(b1,1:2,1));
plot(X',Y','k')

X = node_x(boundary(b2,1:2,1));
Y = node_y(boundary(b2,1:2,1));
plot(X',Y','-r')

%Plotting location of new nodes
scatter(newnode_x,newnode_y,'.k');



end

