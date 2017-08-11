function plot_mesh(mesh,meshfile,field_tmp)
%% plot_mesh.m
%% Author: Timothy Williams
%% Date: 20170811
%% CALL: plot_mesh(mesh,meshfile)
%% mesh = struct eg
%%     Elements: [23616x1 double]
%%      Nodes_x: [4085x1 double]
%%      Nodes_y: [4085x1 double]
%%           id: [4085x1 double]
%% meshfile (optional) is a string
%% - used in plot_coastlines_and_boundaries_c.m if present and exists in path
PLOT_COAST  = 0;
if exist('meshfile','var')
   if exist(meshfile,'file')
      PLOT_COAST  = 1;
   end
end

PLOT_FIELD  = 0;
if exist('field_tmp','var')
   PLOT_FIELD  = 1;
end

sfac  = 1e-3;%convert to km

Nn    = length(mesh.Nodes_x);
Ne    = length(mesh.Elements)/3;
els   = reshape(mesh.Elements,3,Ne)';
xnods = sfac*mesh.Nodes_x(els);
ynods = sfac*mesh.Nodes_y(els);
xels  = mean(xnods,2);
yels  = mean(ynods,2);

%% plot boundaries and coastlines
if PLOT_COAST
   plot_coastlines_and_boundaries_c(meshfile);
   %pause;
   hold on;
end

%% plot boundaries of elements
for j=1:Ne
   xx = [xnods(j,:),xnods(j,1)];
   yy = [ynods(j,:),ynods(j,1)];
   plot(xx,yy,'k');
   hold on;
end

%% plot centres of elements
plot(xels,yels,'.r');
hold off;

%% make full screen to see mesh better...
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

if PLOT_FIELD
   Nd = length(field_tmp);
   if(Nd==Ne)
      % scalar on elements
      v{1}=[field_tmp,field_tmp,field_tmp]';
   elseif(Nd==2*Nn)
      % vector on nodes
      var_mc=field_tmp(mesh.Elements);
      v{1}=reshape(var_mc,[3,Ne]);
      var_mc=field_tmp(mesh.Elements+Nn);
      v{2}=reshape(var_mc,[3,Ne]);
      v{3}=hypot(v{1},v{2});
   elseif(Nd==Nn)
      var_mc=field_tmp(mesh.Elements);
      v{1}=reshape(var_mc,[3,Ne]);
   else
      ss = sprintf('Not the right dimensions - %d - should be %d, %d, or %d',Nd,Ne,Nn,2*Nn);
      error(ss);
   end

   figure;
   patch(xnods',ynods',v{end},'FaceColor','flat','EdgeColor',[0.3 0.3 0.3],'LineWidth',0.2)
   colorbar;
   set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
end
return;
