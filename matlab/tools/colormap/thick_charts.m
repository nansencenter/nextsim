function map = thick_charts(varargin)
%GRAY2RED  colormap used for thick charts
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

map=zeros(300,3); % 1 per cm
map(1:1,:)       = ones(1,1)*w ;
map(2:39,:)      = ones(39-2+1,1)*(0.4*b+0.6*w) ;
map(40:59,:)     = ones(59-40+1,1)*b ;
map(60:79,:)     = ones(79-60+1,1)*(0.5*g+0.5*w);
map(80:99,:)     = ones(99-80+1,1)*g;
map(100:119,:)   = ones(119-100+1,1)*(0.5*y+0.5*w);
map(120:179,:)   = ones(179-120+1,1)*y;
map(180:249,:)   = ones(249-180+1,1)*(r+y)/2;
map(250:299,:)   = ones(299-250+1,1)*r;
map(300,:)       = ones(1,1)*(k+w)/2;
