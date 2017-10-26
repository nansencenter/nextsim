function hh = m_trimesh2D(tri,varargin)
%m_TRIMESH2D Triangular mesh plot in 2D accordingly to the projection
%defined in m_proj.
%m_TRIMESH2D(TRI,LAT,LON) displays the triangles defined in the M-by-3
%face matrix TRI as a mesh.  A row of TRI contains indexes into
%the X,Y, and Z vertex vectors to define a single triangular face.

if (nargin < 1)
     error(message('MATLAB:trimesh:NotEnoughInputs'));
end
d = tri(:,[1 2 3 1])';
lat = varargin{1};
lon = varargin{2};

if nargin == 3
   m_plot(lon(d), lat(d));
else
   m_plot(lon(d), lat(d),varargin{4:end});
end
return;