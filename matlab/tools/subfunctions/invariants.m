function [invar]=invariants(dudx,dudy,dvdx,dvdy)
%function [invar]=invariants(dudx,dudy,dvdx,dvdy)
%compute the invariants div, shear, total deformation rate and vorticity

div=dudx+dvdy;
shear=sqrt( (dudx-dvdy).^2 + (dudy+dvdx).^2 );
eps= sqrt( shear.^2 + div.^2);
vorticity=dvdx-dudy;

invar.div=div;
invar.shear=shear;
invar.eps=eps;
invar.vor=vorticity;

return