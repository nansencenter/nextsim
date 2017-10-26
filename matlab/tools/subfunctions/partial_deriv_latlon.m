function [dudlon,dudlat,dvdlon,dvdlat,scale,area]=partial_deriv_latlon(lon,lat,a,u,v)
%  Authors: Sylvain Bouillon, summer 2013
%
%  GOAL: to compute partial derivatives, using the contour integral method (RGPS user's manual, R. Kwok)
%
%  INPUT:
%  lon, lat: coordinates of the grid
%  u, v: velocity on the grid
%
%  OUTPUT:
%  dudlon, dudlat,...: partial derivatives
%  scale, area : scale (=sqrt(area)) and area of the elements

[nb_lat,nb_lon]=size(u);
dudlon=zeros(nb_lat-1,nb_lon-1);
dudlat=zeros(nb_lat-1,nb_lon-1);
dvdlon=zeros(nb_lat-1,nb_lon-1);
dvdlat=zeros(nb_lat-1,nb_lon-1);
scale=zeros(nb_lat-1,nb_lon-1);
area=zeros(nb_lat-1,nb_lon-1);

for i=1:nb_lat-1
    for j=1:nb_lon-1
        dudlon(i,j)=1/(a*cos(lat(i,j)))*(u(i,j+1)-u(i,j))/(lon(i,j+1)-lon(i,j));
        dvdlon(i,j)=1/(a*cos(lat(i,j)))*(v(i,j+1)-v(i,j))/(lon(i,j+1)-lon(i,j));
        dudlat(i,j)=1/(a*cos(lat(i,j)))*(u(i+1,j)*cos(lat(i+1,j))-u(i,j)*cos(lat(i,j)))/(lat(i+1,j)-lat(i,j));
        dvdlat(i,j)=1/(a*cos(lat(i,j)))*(v(i+1,j)*cos(lat(i+1,j))-v(i,j)*cos(lat(i,j)))/(lat(i+1,j)-lat(i,j));
        area(i,j)=a*cos(lat(i,j))*(lat(i,j+1)-lat(i,j))*(lon(i+1,j)-lon(i,j));
        scale(i,j)=sqrt(area(i,j));
        %         lon_k=[lon(i,j),lon(i+1,j),lon(i+1,j+1),lon(i,j+1),lon(i,j)];
%         lat_k=[lat(i,j),lat(i+1,j),lat(i+1,j+1),lat(i,j+1),lat(i,j)];
%         u_k=[u(i,j),u(i+1,j),u(i+1,j+1),u(i,j+1),u(i,j)];
%         v_k=[v(i,j),v(i+1,j),v(i+1,j+1),v(i,j+1),v(i,j)];
%         tmp_area=0;
%         tmp_dudlon=0;
%         tmp_dudlat=0;
%         tmp_dvdlon=0;
%         tmp_dvdlat=0;
%         for k=1:4
%             tmp_area  =tmp_area+0.5*(lon_k(k)*lat_k(k+1)-lat_k(k)*lon_k(k+1));
%             tmp_dudlon=tmp_dudlon+0.5*(u_k(k+1)+u_k(k))*(lat_k(k+1)-lat_k(k));
%             tmp_dudlat=tmp_dudlat+0.5*(u_k(k+1)+u_k(k))*(lon_k(k+1)-lon_k(k));
%             tmp_dvdlon=tmp_dvdlon+0.5*(v_k(k+1)+v_k(k))*(lat_k(k+1)-lat_k(k));
%             tmp_dvdlat=tmp_dvdlat+0.5*(v_k(k+1)+v_k(k))*(lon_k(k+1)-lon_k(k));
%         end
%         tmp_area=abs(tmp_area);
%         area(i,j)=tmp_area;
%         dudlon(i,j)=tmp_dudlon/tmp_area;
%         dudlat(i,j)=-tmp_dudlat/tmp_area;
%         dvdlon(i,j)=tmp_dvdlon/tmp_area;
%         dvdlat(i,j)=-tmp_dvdlat/tmp_area;
%         scale(i,j)=sqrt(tmp_area);
    end
end