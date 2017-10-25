function h = blue_ice(m)
% BLUE_ICE    Blue-white color map intended for sea-ice concentration maps.
%   It is simply a linear scaling from the blue of PARULA to white.
%   BLUE_ICE(M) returns an M-by-3 matrix containing a "ice" colormap.
%   BLUE_ICE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(ice)
%
%   See also REV_GRIS, ICE

% Einar Olason Tue 21 Jul 07:49:02 CEST 2015

if nargin < 1, m = size(get(gcf,'colormap'),1); end

try
    blue = parula(m);
catch err % We're using an older MATLAB version
    blue = jet(m);
    % blue(1,1) = .2081;
    % blue(1,2) = .1663;
    % blue(1,3) = .5292;
end

r = (0:m-1)'/(m-1) + blue(1,1)*(m-1:-1:0)'/(m-1);
g = (0:m-1)'/(m-1) + blue(1,2)*(m-1:-1:0)'/(m-1);
b = (0:m-1)'/(m-1) + blue(1,3)*(m-1:-1:0)'/(m-1);

h = [r g b];
