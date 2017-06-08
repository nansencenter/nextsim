function [mesh,element]=msh2mat_c(filename,options)
% CALL:  [mesh,element]=msh2mat_c(filename,options)
% Read a .msh file and create mesh and element
% For speed next time, saves mesh and element to a .mat file
% This matfile is not used if options.OVER_WRITE=1 (default=0),
% and not saved if options.SAVE_MATFILE=0 (default=1)
if exist('options','var')
   fields   = fieldnames(options);
   for n=1:length(fields)
      fld   = fields{n};
      eval([fld,'=options.',fld,';']);
   end
   clear options fields fld;
end

%set default options, if not passed in
if ~exist('OVER_WRITE'  ,'var'); OVER_WRITE   = 0; end
if ~exist('SAVE_MATFILE','var'); SAVE_MATFILE = 1; end

matfile  = strrep(which(filename),'.msh','.mat');
if exist(matfile,'file')&~OVER_WRITE
   disp(['Using ',matfile,' instead of ',filename]);
   load(matfile);
   return;
end

fid=fopen(filename);
Stop  = 0;
while ~Stop
   ss    = fgetl(fid);

   %1st check if (eg) the boundary flags are defined
   if strcmp(ss,'$PhysicalNames')
      Nflags   = fscanf(fid, '%d\n',1);
      for n=1:Nflags
         tc = strsplit(fgetl(fid));
         mesh.flags.(tc{3}(2:end-1))   = str2num(tc{2});
      end
      ss    = fgetl(fid);
      ss    = fgetl(fid);
   end
   Stop  = strcmp(ss,'$Nodes');
end
%fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid)
%fgetl(fid), fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid);

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
Stop  = 0;
while ~Stop
   ss    = fgetl(fid);
   Stop  = strcmp(ss,'$Elements');
end
%fgetl(fid); fgetl(fid); fgetl(fid);
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
fclose(fid);

% reduction of the size of element.num_node and mesh.boundary.from_msh
element.num_node      =element.num_node      (1:ind_el-1,:);
mesh.boundary.from_msh=mesh.boundary.from_msh(1:ind_bc-1,:);

% number of elements
mesh.Ne=ind_el-1;
if SAVE_MATFILE
   disp(['Saving ',matfile]);
   save(matfile,'mesh','element');
end
return

