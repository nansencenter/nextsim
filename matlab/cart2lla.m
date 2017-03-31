% CART2LLA - convert cartesian coordinates given on a sphere to latitude, longitude
%
% USAGE:
% [lat,lon,radius] = cart2lla(x,y,z)
%
% lat = latitude (radians)
% lon = longitude (radians)
% radius (m)
% x = X-coordinate on a sphere (m)
% y = Y-coordinate on a sphere (m)
% z = Z-coordinate on a sphere (m)
%
% Notes: (1) This function assumes that coordinates are given on a sphere.
%        (2) Inputs may be scalars, vectors, or matrices of the same
%            size and shape. Outputs will have that same size and shape.
%        (3) Tested but no warranty; use at your own risk.
%        (5) Sylvain Bouillon, August 2013

function [lat,lon,radius] = cart2lla(x,y,z)

% calculations:
radius = sqrt(x.^2+y.^2+z.^2);
lat = asin(z./radius);
lon = atan2(y,x);

% return lon in range [0,2*pi)
lon = mod(lon,2*pi);

return
