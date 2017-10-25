function map = metno_charts(varargin)
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

map=zeros(100,3);
map(1:1,:)      = ones(1,1)*w ;
map(2:10,:)     = ones(10-2+1,1)*(0.4*b+0.6*w) ;
map(11:40,:)    = ones(40-11+1,1)*(0.5*g+0.5*w);
map(41:70,:)    = ones(70-41+1,1)*y;
map(71:90,:)    = ones(90-71+1,1)*(r+y)/2;
map(91:100,:)   = ones(100-91+1,1)*r;
%map(91:99,:)   = ones(99-91+1,1)*r;
%map(100,:)      = ones(1,1)*(k+w)/2;
