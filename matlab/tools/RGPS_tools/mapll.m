function [x,y] = mapll(lon,lat,sgn,delta,slat,wgs84)
%    
% NAME:       mapll
% 
% PURPOSE:    This subroutine converts from geodetic latitude and
%             longitude to Polar Stereographic (X,Y) coordinates for
%             the polar regions.
%     
% CATEGORY:   Misc
%     
% USAGE:      [x,y] = mapll(lon,lat[,Options])
% 
% INPUTS:     lon (Float)        = Geodetic Longitude (degrees, 0 to 360)
%             lat (Float)        = Geodetic Latitude  (degrees, +90 to -90)
%             sgn                = Sign of latitude (automatic detected
%                                  if not given)
%                                  +1 : north latitude
%                                  -1 : south latitude
%             delta (Float)      = Angle of rotation of the grid to
%                                  Greenwich 0 degree longitude.
%                                  hemisphere. [default: NSIDC grid].
%             slat  (Float)      = Latitude for the grid with no distortion
%                                  (default: 70.), never set to zero.
% 			        (STRING)= 'zero' for slat = 0.
%             wgs84              = 1: WGS84 ellipsoid; 0: Hough ellipsoid
% 
% 
% OUTPUTS:         x   (Float) = Polar Stereographic X Coordinate (km)
%                  y   (Float) = Polar Stereographic Y Coordinate (km)
%     
% 
% AUTHOR:          The original equations are from Snyder, J. P.,
%                  1982,  Map Projections Used by the U.S. Geological Survey,
%                  Geological Survey Bulletin 153llssmi1.m2, U.S. Government Printing 
%                  Office.  See JPL Technical Memorandum 3349-85-101 for 
%                  further details.
% 
%                  Original FORTRAN program written by C. S. Morris,
%                  April 1985, Jet Propulsion Laboratory, California
%                  Institute of Technology
% 
%                  IDL conversion by Helmut Schottmueller, September 1995,
%                  Institute for environmental physics, University of Bremen
% 
%                  Error correction for RHO: 1997-11-03 by LML
%                  Added some options: 2003-10-06 Gunnar Spreen
%                  More options:       2005-04-08 Gunnar Spreen
%                  automatic detect sgn 2005-06-16 Gunnar Spreen
%                  allow delta to be 0  Nov 2006 Gunnar Spreen
%
%                  Matlab conversion   Mar 2009; Gunnar Spreen
% -

% if not set 'sgn' use largest lat for sgn
if nargin < 3
    switch 1;
        case max(max(lat)) >= 0;
            sgn = 1;
        case max(max(lat)) < 0;
            sgn = -1;
        otherwise;
            disp('Can not detect hemisphere (sgn), please specify!');
            return;
    end;
end;

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

% slat must be positiv, correct this if not the case
slat = abs(slat);

londelta = lon-delta;

% Compute X and Y in grid coordinates.
t = tand(45 - lat./2)./((1-e.*sind(lat))./(1+e.*sind(lat))).^(e/2);

if (abs(90-slat) < 1e-5)
    rho = 2.*re.*t./sqrt((1+e)^(1+e)*(1-e)^(1-e));
else
    tc  = tand(45-(slat/2))/((1-e*sind(slat))/(1+e*sind(slat)))^(e/2);
    mc  = cosd(slat)./sqrt(1-e2*(sind(slat)^2));
    rho = re*mc.*t./tc;
end
x = rho.*sgn.*sind(sgn*londelta);
y = -rho.*sgn.*cosd(sgn*londelta);

ind = find(lat >= 90);
if(~isempty(ind))
   x(ind)=0;
   y(ind)=0;
end



