function map = gray2red(varargin)
%GRAY2RED  colormap showing a gradient from gray to red.
%
%   MAP = gray2red
%   return a 256x3 array which can be used with colormap.
%
%   Example
%   colormap(gray2red);
%
%   See also
%   colormap
%
%

r = zeros(256,1);

r = zeros(256,1);
r(1:256) = linspace(100, 75, 256);

g = zeros(256,1);
g(1:256) = linspace(0, 75, 256);

b = zeros(256,1);
b(1:256) = linspace(0,75,256);

map = [flipud(r) flipud(g) flipud(b)]/100;
