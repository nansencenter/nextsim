function map = blue2red(varargin)
%BLUE2RED  colormap showing a gradient from blue to green to red.
%
%   MAP = blue2red
%   return a 256x3 array which can be used with colormap.
%   map(1,:)    corresponds to blue color,
%   map(64,:)   corresponds to cyan color,
%   map(128,:)  corresponds to green color,
%   map(192,:)  corresponds to yellow color,
%   map(256,:)  corresponds to red color,
%   and all indices inbetween are gradient between the 2 extremes colors.
%
%
%   Example
%   colormap(blue2red);
%
%   See also
%   colormap
%
%
% ------
% Author: David Legland
% e-mail: david.legland@jouy.inra.fr
% Created: 2006-05-24
% Copyright 2006 INRA - CEPIA Nantes - MIAJ (Jouy-en-Josas).


r = zeros(256,1);
r(1:128) = linspace(100, 75, 128);
r(129:256)=linspace(75, 0, 128);

g = zeros(256,1);
g(1:128) = linspace(0, 75, 128);
% g(1:128) = linspace(0, 100, 128);
g(129:256)=linspace(75, 0, 128);

b = zeros(256,1);
b(1:128) = linspace(0,75,128);
b(129:256) = linspace(75,100,128);

% r = zeros(256,1);
% r(1:106) = linspace(58, 100, 106);
% r(107:146)=100;
% r(147:195)=linspace(100, 0, 49);
% 
% g = zeros(256,1);
% g(61:106) = linspace(20, 100, 46);
% g(107:146) = 100;
% g(147:256)=linspace(100, 16, 110);
% 
% b = zeros(256,1);
% b(106:147)=3;
% b(148:256) = linspace(71,100,109);

map = [flipud(r) flipud(g) flipud(b)]/100;
