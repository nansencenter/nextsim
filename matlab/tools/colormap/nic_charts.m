function map = nic_charts(varargin)
%GRAY2RED  colormap used for NIC ice charts
%
%   MAP = NIC_charts
%   return a 101x3 array which can be used with colormap.
%
%   Example
%   colormap(NIC_charts);
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

map=zeros(101,3);
map(1:18,:) = ones(18,1)*b ;
map(19:80,:) = ones(80-19+1,1)*y;
map(81:101,:) = ones(101-81+1,1)*r;
