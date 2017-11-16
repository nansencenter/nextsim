function h = ice(m)
%ICE    Blue-white color map intended for se-ice concentration maps. It is
%   based on the HOT colormap (like REV_GRIS), but goes from a dark blue color to
%   white, instead of black-blue-white.
%   ICE(M) returns an M-by-3 matrix containing a "ice" colormap.
%   ICE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(ice)
%
%   See also REV_GRIS, BLUE_ICE

% Einar Olason Tue 21 Jul 07:49:02 CEST 2015

if nargin < 1, m = size(get(gcf,'colormap'),1); end
n = fix(3/8*m);

b = [min(1, 0.5+(1:2*n)'/(2*n)); ones(m-2*n,1)];
g = [(1:2*n)'/(2*n); ones(m-2*n,1)];
r = [zeros(2*n,1); (1:m-2*n)'/(m-2*n)];

h = [r g b];
