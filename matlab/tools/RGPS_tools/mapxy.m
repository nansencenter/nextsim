function [lon,lat] = mapxy(x,y,sgn,delta,slat,wgs84)
% +  
% NAME:       mapxy
%    
% PURPOSE:    This subroutine converts from Polar Stereographic (X, Y) coordinates
%             to geodetic latitude and longitude for the polar regions.
%     
% USAGE:      [lon,lot] = mapxy(x,y [,options])
% 
% INPUTS:     x          (Float) = Polar Stereographic X Coordinate (km)
%             y          (Float) = Polar Stereographic Y Coordinate (km)
%             sgn                = Sign of latitude
%                                 +1 : north latitude [default]
%                                 -1 : south latitude
%            delta      (Float) = Angle of rotation of the grid
%                                 to Greenwich 0 degree longitude
%                                 [default: NSIDC grid]
%             slat       (Float) = Latitude of true scale for the grid (default: 70.).
%                                  Never set zero.
%            wgs84              = 1: use WGS84 ellipsoid; 0: Hough ellipsoid
% 
% 
% OUTPUTS:         lon  (Float) = Geodetic Longitude
%                  lat  (Float) = Geodetic Latitude
%     
% 
% AUTHOR:          The original equations are from Snyder, J. P.,
%                  1982,  Map Projections Used by the U.S. Geological Survey,
%                  Geological Survey Bulletin 1532, U.S. Government Printing 
%                  Office.  See JPL Technical Memorandum 3349-85-101 for 
%                  further details.
% 
%                  Original FORTRAN program written by C. S. Morris,
%                  April 1985, Jet Propulsion Laboratory, California
%                  Institute of Technology
% 
%                  IDL conversion by Helmut Schottmueller, September 1995,
%                  Institute for environmental physics, University of
%                  Bremen
% 
%                  Added some options: 2003-10-06 Gunnar Spreen
%                  WGS84:              2005-09-20 Gunnar Spreen
%                  Matlab conversion   Mar 2009; Gunnar Spreen
% -

% x and y must be at least floats
if ~(isfloat(x) && isfloat(y))
    y=double(y);
    x=double(x);
end

% if not set 'sgn' use northern hemisphere
if nargin < 3, sgn = 1; end

if nargin < 4
    switch sgn
        case 1, delta = -45.;
        case -1, delta = 0;
    end
end

% Standard latitude of true scale for the SSM/I grid
if nargin < 5, slat = 70.; end

if (nargin < 6 || wgs84 == 0)   % default is the Hough ellipsoid
    % Radius of the earth in kilometers (A)
    re = 6378.273;
    % Eccentricity of the Hughes ellipsoid squared
    e2 = 0.006693883;       % => B=6356.8894
    % Eccentricity of the Hughes ellipsoid
    e  = sqrt(e2);
else
    if wgs84    % use WGS84 ellipsoid
        re = 6378.137;
        e  = sqrt(1-(6356.7523^2/re^2));    % => B=6356.7523
        e2 = e^2;
    else
        disp('No ellipsoid defined!')
        return
    end
end

% slat must be positiv, correct this
slat = abs(slat);

sl = slat * pi/180;

rho = sqrt(x.^2 + y.^2);
cm = cos(sl) / sqrt(1 - e2 * (sin(sl)^2));
tc = tan(pi/4 - sl/2) / ((1 - e*sin(sl))/(1 + e*sin(sl)))^(e/2);
if (abs(slat-90) < 1e-5)
    t = rho.* sqrt((1+e)^(1+e).* (1-e)^(1-e))./2./re;
else
    t = rho.* tc./(re*cm);
end
chi = pi/2 - 2*atan(t);
lat = chi + (e2/2 + 5* e2^2 /24 + e2^3 /12) * ...
      sin(2*chi) + ((7* e2^2 /48) + (29* e2^3 /240)) * ...
      sin(4*chi) + (7* e2^3 /120) * sin(6*chi);
lat = sgn.*lat;
lon = atan2(sgn.*x, -sgn.*y);
lon = sgn.*lon;

lon = lon.* 180/pi;
lat = lat.* 180/pi;
lon = lon + delta;

ind = find(rho <= 0.001);
if ~isempty(ind)
    lat(ind) = 90 * sgn;
    lon(ind) = 0;
end