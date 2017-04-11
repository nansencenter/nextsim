function diff_final(outdir1,outdir2,step,PLOT_ANOMALIES)

if ~exist('step','var')
   step  = 1000;%final results
end
if ~exist('PLOT_ANOMALIES','var')
   PLOT_ANOMALIES  = 0;%don't plot anomalies
end

try
    [mesh1, data1] = neXtSIM_bin_revert(outdir1,[],step);
    [mesh2, data2] = neXtSIM_bin_revert(outdir2,[],step);
catch
    error('Cannot read data');
end

[Mesh1,Data1] = convert_mesh(mesh1, data1);
[Mesh2,Data2] = convert_mesh(mesh2, data2);


% check if on the same mesh
SAME_MESH   = 0;
if Mesh1.Nn==Mesh2.Nn & Mesh1.Ne==Mesh2.Ne
   % meshes could be the same
   dist  = hypot(Mesh1.Nodes_x-Mesh2.Nodes_x,Mesh1.Nodes_y-Mesh2.Nodes_y);%distance in m
   if max(dist)<1e-3
      SAME_MESH   = 1;
   end
end


if SAME_MESH
   disp('Results are on the same mesh');
   %
   disp(' ');
   disp('Comparing scalars...');
   compare_fields(Data1.Scalars,Data2.Scalars);
   %
   disp(' ');
   disp('Comparing elements...');
   compare_fields(Data1.Elements,Data2.Elements);
   %
   disp(' ');
   disp('Comparing nodes...');
   compare_fields(Data1.Nodes,Data2.Nodes);
else
   disp('Results are on different meshes');
   disp(' - interpolating onto a grid');
   Grid  = make_grid(Mesh1,Mesh2);
   save Grid Grid
   %
   disp(' ');
   disp('Comparing scalars...');
   %make a high-res grid to interpolate onto
   anomaly.Scalars   = compare_fields(Data1.Scalars,Data2.Scalars);

   %interp and compare elements
   disp(' ');
   disp('Comparing elements...');
   USE_ISSM = 1;
   if USE_ISSM
      Gdata1   = interp_ISSM(Grid,Mesh1,Data1.Elements);
      Gdata2   = interp_ISSM(Grid,Mesh2,Data2.Elements);
   else
      Gdata1   = interp_griddata(Grid,Mesh1.Elements_x,Mesh1.Elements_x,Data1.Elements);
      Gdata2   = interp_griddata(Grid,Mesh2.Elements_x,Mesh2.Elements_x,Data2.Elements);
   end
   anomaly.Elements  = compare_fields(Gdata1,Gdata2);

   %interp and compare nodes
   disp(' ');
   disp('Comparing nodes...');
   if USE_ISSM
      Gdata1   = interp_ISSM(Grid,Mesh1,Data1.Nodes);
      Gdata2   = interp_ISSM(Grid,Mesh2,Data2.Nodes);
   else
      Gdata1   = interp_griddata(Grid,Mesh1.Nodes_x,Mesh1.Nodes_x,Data1.Nodes);
      Gdata2   = interp_griddata(Grid,Mesh2.Nodes_x,Mesh2.Nodes_x,Data2.Nodes);
   end
   anomaly.Nodes  = compare_fields(Gdata1,Gdata2);
end

if PLOT_ANOMALIES
   %plot anomalies
   if ~SAME_MESH
      plot_fields_grid(anomaly.Elements);
      plot_fields_grid(anomaly.Nodes);
   else
      error('Plotting of anomalies on mesh not implemented');
   end
end


return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_fields_grid(Data,Grid)

IMSHOW   = (~exist('Grid','var'));
fields   = fieldnames(Data);
if IMSHOW
   for j=1:length(fields)
      figure;
      fld   = fields{j};
      data  = flipud(Data.(fld)');
      imshow(data);
      title(strrep(fld,'_','-'));
      colorbar;

      caxis([min(data(:)),max(data(:))]);
      colormap('jet');
   end
else
   error('Plotting using grid info not implemented');
   %TODO plot with pcolor and then use plot_coastlines_and_boundaries_c.m
   %TODO OR plot with m_proj
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function anomaly  = compare_fields(data1,data2)
err = 0;
fields=fieldnames(data1);
for i=1:length(fields)
   fld   = fields{i};
   anomaly.(fld) = abs(data2.(fld) - data1.(fld));
   %%
   jok   = find(~isnan(anomaly.(fld)));
   err1  = max(abs(anomaly.(fld)(jok)));
   disp(['Maximum absolute error in ',fld,' = ',num2str(err1,'%f')])
   disp(['  Minimums: ',num2str(min(data1.(fld)(:))),' | ',num2str(min(data2.(fld)(:)))])
   disp(['  Maximums: ',num2str(max(data1.(fld)(:))),' | ',num2str(max(data2.(fld)(:)))])
   %err   = max(err,err1);
end

%disp(['Maximum error: ' num2str(err)])
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mesh,Data] = convert_mesh(mesh,data)
% input:
% mesh:
%   Elements: [126348x1 double]
%    Nodes_x: [21912x1 double]
%    Nodes_y: [21912x1 double]
%         id: [21912x1 double]
%
% output:
% Mesh = 
%                   Elements: [126348x1 double]
%                    Nodes_x: [21912x1 double]
%                    Nodes_y: [21912x1 double]
%                         id: [21912x1 double]
%               Element_area: [42116x1 double]  %element area as calculated by nextsim
%                AllMinAngle: [42116x1 double]
%          PreviousNumbering: [21912x1 double]
%          M_dirichlet_flags: [1349x1 double]
%                         Nn: 21912
%                         Ne: 42116
%     Element_area_projected: [42116x1 double]  %calculate element area in a Euclidean way from (x,y) coordinates
%                 Elements_x: [42116x1 double]
%                 Elements_y: [42116x1 double]
%                 resolution: 8.5772e+03
%                       xmin: 1.4200e+05
%                       xmax: 1.5180e+06
%                       ymin: -1.8109e+06
%                       ymax: -1.7086e+05

% coords of nodes
Mesh  = mesh;
clear mesh;

% get some mesh stuff from data
flds  = {'Element_area','AllMinAngle','PreviousNumbering','M_dirichlet_flags'};
for j=1:length(flds)
   fld   = flds{j};
   if isfield(data,fld)
      Mesh.(fld)  = data.(fld);
      data  = rmfield(data,fld);
   end
end

% size of mesh
Mesh.Nn  = length(Mesh.Nodes_x);
Mesh.Ne  = length(Mesh.Elements)/3;

% get the 3 nodes surrounding each element
var_mx   = reshape(Mesh.Nodes_x(Mesh.Elements),[3,Mesh.Ne]);
var_my   = reshape(Mesh.Nodes_y(Mesh.Elements),[3,Mesh.Ne]);

% coords of elements
Mesh.Element_area_projected   = polyarea(var_mx,var_my)';
Mesh.Elements_x               = mean(var_mx)';
Mesh.Elements_y               = mean(var_my)';

% resolution from area = b*h/2 = h^2/2/sqrt(3)
% for an equilateral triangle (a,b,c = 1,sqrt(3),2)
hsq               = 2*sqrt(3)*min(Mesh.Element_area);%h^2
Mesh.resolution   = sqrt(hsq);

%other info
Mesh.xmin   = min(Mesh.Nodes_x);
Mesh.xmax   = max(Mesh.Nodes_x);
Mesh.ymin   = min(Mesh.Nodes_y);
Mesh.ymax   = max(Mesh.Nodes_y);

Data  = sort_data(Mesh,data);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = make_grid(mesh1,mesh2)
% make grid to cover area spanned by both meshes
% inputs:
% mesh1,mesh2 = 
%                   Elements: [126348x1 double]
%                    Nodes_x: [21912x1 double]
%                    Nodes_y: [21912x1 double]
%                         id: [21912x1 double]
%               Element_area: [42116x1 double]  %element area as calculated by nextsim
%                AllMinAngle: [42116x1 double]
%          PreviousNumbering: [21912x1 double]
%          M_dirichlet_flags: [1349x1 double]
%                         Nn: 21912
%                         Ne: 42116
%     Element_area_projected: [42116x1 double]  %calculate element area in a Euclidean way from (x,y) coordinates
%                 Elements_x: [42116x1 double]
%                 Elements_y: [42116x1 double]
%                 resolution: 8.5772e+03
%                       xmin: 1.4200e+05
%                       xmax: 1.5180e+06
%                       ymin: -1.8109e+06
%                       ymax: -1.7086e+05
% output:
% Grid = 
%     xmin: 1.4200e+05
%     ymin: -1.8109e+06
%     xmax: 1.5180e+06
%     ymax: -1.7086e+05
%       nx: 325
%       dx: 4.2338e+03
%       ny: 387
%       dy: 4.2377e+03
%        x: [325x1 double]
%        y: [387x1 double]
%        Y: [325x387 double]
%        X: [325x387 double]
% x increases down rows, y increases across columns

G.xmin  = min(mesh1.xmin,mesh2.xmin);
G.ymin  = min(mesh1.ymin,mesh2.ymin);
G.xmax  = max(mesh1.xmax,mesh2.xmax);
G.ymax  = max(mesh1.ymax,mesh2.ymax);

res   = .5*min(mesh1.resolution,mesh2.resolution);
G.nx    = ceil((G.xmax-G.xmin)/res);
G.dx    = (G.xmax-G.xmin)/G.nx;
G.ny    = ceil((G.ymax-G.ymin)/res);
G.dy    = (G.ymax-G.ymin)/G.ny;

G.x         = G.xmin+.5*G.dx+(1:G.nx)'*G.dx;
G.y         = G.ymin+.5*G.dy+(1:G.ny)'*G.dy;
[G.Y,G.X]   = meshgrid(G.y,G.x);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Data  = sort_data(mesh,data)

fields=fieldnames(data);
for i=1:length(fields)
    fld  = fields{i};
    dat  = data.(fld);
    if length(dat)==1
       % global scalar
       Data.Scalars.(fld)  = dat;
    elseif length(dat)==mesh.Ne
       % scalar on elements
       Data.Elements.(fld) = dat;
    elseif length(dat)==mesh.Nn
       % scalar on nodes
       Data.Nodes.(fld) = dat;
    elseif length(dat)==2*mesh.Nn
       % vector on nodes
       Data.Nodes.([fld,'_x'])   = dat(1:mesh.Nn);
       Data.Nodes.([fld,'_y'])   = dat(mesh.Nn+1:2*mesh.Nn);
    else
       {fld,dat,mesh.Ne,mesh.Nn}
       %Data.Other.(fld) = dat;
       error('Not the right dimensions')
    end
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grid_data = interp_griddata(gridprams,x_mesh,y_mesh,mesh_data,method)
%% method:
%%   'nearest'   - Nearest neighbor interpolation
%%    'linear'    - Linear interpolation (default)
%%    'natural'   - Natural neighbor interpolation
%%    'cubic'     - Cubic interpolation (2D only)

if ~exist('method','var')
   method   = 'linear';
end
fields      = fieldnames(mesh_data);
nvar        = length(fields);
grid_data   = zeros(gridprams.nx,gridprams.ny,nvar);
xx          = gridprams.X(:);
yy          = gridprams.Y(:);
for j=1:nvar
   fld   = fields{j}
   {xx,yy,gridprams.nx*gridprams.ny}
   data              = griddata(x_mesh,y_mesh,mesh_data.(fields{j}),xx,yy,method);
   {data}
   grid_data.(fld)   = reshape(data,gridprams.nx,gridprams.ny);
end

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grid_data = interp_ISSM(gridprams,mesh,mesh_data,method)

xmin     = min(gridprams.x);  % WIM grid xmin
ymax     = max(gridprams.y);  % WIM grid ymax
xposting = gridprams.dx;      % res in x dirn (km)
yposting = gridprams.dy;      % res in y dirn (km)
nlines   = gridprams.ny;      % no of cols in WIM grid
ncols    = gridprams.nx;      % no of rows in WIM grid

inputs.x             = mesh.Nodes_x;
inputs.y             = mesh.Nodes_y;
inputs.index         = reshape(mesh.Elements,[3,mesh.Ne])';
inputs.xmin          = xmin;
inputs.ymax          = ymax;
inputs.xposting      = xposting;
inputs.yposting      = yposting;
inputs.nlines        = double(nlines); %% inputs need to be doubles
inputs.ncols         = double(ncols);  %% inputs need to be doubles
inputs.default_value = NaN;            %% value if outside mesh

fields   = fieldnames(mesh_data);
nvar     = length(fields);
for j=1:nvar
   fld         = fields{j};
   inputs.data = mesh_data.(fld);

   [x_m,y_m,Gdata]   =...
      InterpFromMeshToGrid(inputs.index,inputs.x,inputs.y,inputs.data,...
                           inputs.xmin,inputs.ymax,inputs.xposting,inputs.yposting,...
                           inputs.nlines,inputs.ncols,...
                           inputs.default_value);

   grid_data.(fld)   = Gdata';

   if 0
      % Do a test plot as a check
      {min(Gdata(:)),max(Gdata(:))}
      imshow(flipud(grid_data.(fld)'))
      colorbar;
      error('hey');
   end

   if 0
      {grid_data.(fld),gridprams.nx,gridprams.ny},pause
      tst_xmin = [xmin,min(gridprams.X(:)),min(x_m)]
      tst_xmax = [max(gridprams.x),max(x_m)]
      tst_ymin = [min(gridprams.y),min(y_m)]
      tst_ymax = [ymax,max(gridprams.y),max(y_m)]
      tst_dx   = [gridprams.dx,x_m(2)-x_m(1)]
      tst_dy   = [gridprams.dy,y_m(2)-y_m(1)]
      {gridprams.x,x_m}
      {gridprams.y,y_m}
      error('hey')
   end
end
