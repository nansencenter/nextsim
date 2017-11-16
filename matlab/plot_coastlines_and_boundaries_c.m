function plot_coastlines_and_boundaries_c(mesh_filename,RGPS_projection)
switch nargin
    case 1
        RGPS_projection=false;
end

% The files defining the domain are on /Data/sim/data/mesh and have to be
% in your matlab path

% fprintf('Coastlines are being loaded from the file ''%s'' ...\r',mesh_filename)
if(strcmp(mesh_filename(end-3:end),'.msh'))
    mesh=msh2mat_c(mesh_filename);
elseif(strcmp(mesh_filename(end-3:end),'.mat'))
    load(mesh_filename);
end;

if isfield(mesh,'flags')
   flag_boundary_fix=mesh.flags.coast;
   flag_boundary_free=mesh.flags.open;
   % mesh.flags comes from running msh2mat_c.m on .msh file,
   % which reads the following lines (with the correct values for open and coast)
   % underneath $EndMeshFormat:
   %
   % $PhysicalNames
   % 2
   % 1 10001 "open"
   % 1 10000 "coast"
   % $EndPhysicalNames
   %
   % If they are not present, could manually edit the .msh file, then rerun
else
   flags   = sort(unique(mesh.boundary.from_msh(:,3)));

   %may need to be changed depending on the meshfile (most have free = fix+1)
   flag_boundary_fix=flags(1);
   flag_boundary_free=flags(2);
end

boundary   = mesh.boundary.from_msh;

%Selecting closed boundaries
fix = find(flag_boundary_fix==boundary(:,3));

%Selecting free boundaries
free   = [];
for loop_i=1:length(flag_boundary_free)
    fbf = flag_boundary_free(loop_i);
    free= [free;find(fbf==boundary(:,3))];
end

sfac     = 1e-3;%km
R_earth  = 6378.273;%radius of earth (km)
if 1
   %use lon/lat
   node_lat   = mesh.node.lat;
   node_lon   = mesh.node.lon;
   closed_boundaryLat  = node_lat(boundary(fix ,1:2,1))';
   closed_boundaryLon  = node_lon(boundary(fix ,1:2,1))';
   free_boundaryLat = node_lat(boundary(free ,1:2,1))';
   free_boundaryLon = node_lon(boundary(free ,1:2,1))';

   mppfile = which('NpsNextsim.mpp');
   if strcmp(mppfile,'')
      error('add path to mppfile');
   end


<<<<<<< HEAD
   % projection
   if(~RGPS_projection)
       [closed_boundaryX,closed_boundaryY]= mapx_forward(mppfile,closed_boundaryLon(:)',closed_boundaryLat(:)');
       [free_boundaryX,free_boundaryY]= mapx_forward(mppfile,free_boundaryLon(:)',free_boundaryLat(:)');
       
       %reshape and scaling to km
       closed_boundaryX = reshape(sfac*closed_boundaryX,size(closed_boundaryLon));
       closed_boundaryY = reshape(sfac*closed_boundaryY,size(closed_boundaryLon));
       free_boundaryX   = reshape(sfac*free_boundaryX  ,size(free_boundaryLon));
       free_boundaryY   = reshape(sfac*free_boundaryY  ,size(free_boundaryLon));
   else
   % projection used when reading the RGPS Lagrangian data    
       m_proj('Stereographic','lon',-45,'lat',90,'radius',60);
       
       [closed_boundaryX,closed_boundaryY]=m_ll2xy(closed_boundaryLon,closed_boundaryLat);
       closed_boundaryX=double(closed_boundaryX)*R_earth;
       closed_boundaryY=double(closed_boundaryY)*R_earth;
       
       [free_boundaryX,free_boundaryY]=m_ll2xy(free_boundaryLon,free_boundaryLat);
       free_boundaryX=double(free_boundaryX)*R_earth;
       free_boundaryY=double(free_boundaryY)*R_earth;
   end
   
else
   %use x/y
   node_x   = sfac*mesh.node.x;
   node_y   = sfac*mesh.node.y;
   closed_boundaryX  = node_x(boundary(fix ,1:2,1))';
   closed_boundaryY  = node_y(boundary(fix ,1:2,1))';
   free_boundaryX    = node_x(boundary(free ,1:2,1))';
   free_boundaryY    = node_y(boundary(free ,1:2,1))';
end

%Plotting closed mesh boundaries (e.g. coastlines)
plot(closed_boundaryX,closed_boundaryY,'Color',[0.3 0.3 0.3],'LineWidth',0.2);
hold on;

%Plotting open mesh boundaries
plot(free_boundaryX,free_boundaryY,'g','LineWidth',2);
hold off;

return;
