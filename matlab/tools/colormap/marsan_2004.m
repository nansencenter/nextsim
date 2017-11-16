function map = marsan_2004(varargin)
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

y=[1 1 0];
m=[1 0 1];
c=[0 1 1];
r=[1 0 0];
g=[0 1 0];
b=[0 0 1];
w=[1 1 1];
k=[0 0 0];

map = [y; c; b; m];
