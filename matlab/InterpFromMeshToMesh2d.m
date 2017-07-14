function nterp_out = InterpFromMeshToMesh2d(index_data, x_data, y_data, data, x_interp, y_interp, isdefault, defaultvalue)

% INTERFROMMESHTOMESH2D - interpolation from a 2d triangular mesh onto a list of point
% 
%    This function is a multi-threaded mex file that interpolates a field
%    defined on a Delaunay triangulation onto a list of point
% 
%    Usage:
%          data_interp=InterpFromMeshToMesh2d(index,x,y,data,x_interp,y_interp);
%       or data_interp=InterpFromMeshToMesh2d(index,x,y,data,x_interp,y_interp,OPTIONS);
% 
%       index             : index of the mesh where data is defined
%       x,y               : coordinates of the nodes where data is defined
%       data              : matrix holding the data to be interpolated onto the mesh. (one column per field)
%       x_interp,y_interp : coordinates of the points onto which we interpolate.
%       data_interp       : vector of mesh interpolated data.
%       Available options :
%          - 'default' : default value if point is outsite of triangulation (instead of linear interolation)
% 
%    Example:
%       load('temperature.mat');
%       md.initialization.temperature=InterpFromMeshToMesh2d(index,x,y,temperature,md.mesh.x,md.mesh.y);
%       md.initialization.temperature=InterpFromMeshToMesh2d(index,x,y,temperature,md.mesh.x,md.mesh.y,'default',253);
% 
% 
