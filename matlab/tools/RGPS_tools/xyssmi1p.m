function [lon,lat]=xyssmi1p(x,y,clon)
%------------------------------------------------------
% function [lon,lat]=xyssmi1p(x,y,clon)
%
% function to convert [X,Y] in "secant stereographic projection"
% to [lon,lat].  The secant plane intersects at 70N and centered
% at [90N,-45E] or [90N,clon].
% Program taken from mapx package at NSIDC
% 
% INPUT: X,Y in meters		[Mx1] and [Mx1]
% OUTPUT: lon,lat in degress	[Mx1] and [Mx1]
%
% ATN 04.17.07
%------------------------------------------------------
crd=pi/180;			% conversion rad -> deg
cdr=180/pi;
R=6378.273;			% Radius of the Earth (km)
e2=0.006693883;			% Eccentricity square
e=sqrt(e2);
e4=e^4;
e6=e^6;
e8=e^8;
slat=70.0;			% Standard Parallel (degrees)

lam=atan2(x,-y);

rhoi=sqrt(x.^2+y.^2)./1e3;	% in km to be consisten with R unit
numi=rhoi.*tan(pi/4-slat*crd/2)./((1-e*sin(slat*crd))./(1+e*sin(slat*crd))).^(e/2) ;
demi=R*cos(slat*crd)/sqrt(1-e2*(sin(slat*crd)).^2);
ti=numi./demi;
chi=pi/2-2.*atan(ti);
sin2chi=sin(2*chi);
sin4chi=sin(4*chi);
sin6chi=sin(6*chi);

phi=chi+ sin2chi.*e2./2    + sin2chi.*5.*e4./24 ...
       + sin2chi.*e6./12   + sin2chi.*13.*e8./360 ...
       + sin4chi.*7.*e4/48 + sin4chi.*29.*e6./240 ...
       + sin4chi.*811.*e8./11520 ...
       + sin6chi.*7.*e6./120 ...
       + sin6chi.*81.*e8./1120 ...
       + sin(8*chi).*4279.*e8./161280;

if(slat>0);
  lon=lam.*cdr+clon;
  lat=phi.*cdr;
else;
  lon=-lam.*cdr+clon;
  lat=-phi.*cdr;
end;

return
