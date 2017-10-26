function [dudx,dudy,dvdx,dvdy,scale,area,tri_regular]=partial_deriv_regular(x,y,u,v)
%  function [dudx,dudy,dvdx,dvdy,scale,area]=partial_deriv_regular(xy,uv)
%  Authors: Sylvain Bouilx, summer 2013
%
%  GOAL: to compute partial derivatives, using the contour integral method (RGPS user's manual, R. Kwok)
%
%  INPUT:
%  x, y: coordinates of the grid
%  u, v: velocity on the grid
%
%  OUTPUT:
%  dudx, dudy,...: partial derivatives
%  scale, area : scale (=sqrt(area)) and area of the elements

[nb_x,nb_y]=size(u);

dudx=zeros(nb_x-1,nb_y-1);
dudy=zeros(nb_x-1,nb_y-1);
dvdx=zeros(nb_x-1,nb_y-1);
dvdy=zeros(nb_x-1,nb_y-1);
scale=zeros(nb_x-1,nb_y-1);
area=zeros(nb_x-1,nb_y-1);

for i=1:nb_x-1
    for j=1:nb_y-1
%         dudx(i,j)=(u(i+1,j)-u(i,j))/(x(i+1,j)-x(i,j));
%         dvdx(i,j)=(v(i+1,j)-v(i,j))/(x(i+1,j)-x(i,j));
%         dudy(i,j)=(u(i,j+1)-u(i,j))/(y(i,j+1)-y(i,j));
%         dvdy(i,j)=(v(i,j+1)-v(i,j))/(y(i,j+1)-y(i,j));
%         area(i,j)=(y(i,j+1)-y(i,j))*(x(i+1,j)-x(i,j));
%         scale(i,j)=sqrt(area(i,j));
        x_k=[x(i,j),x(i+1,j),x(i+1,j+1),x(i,j+1),x(i,j)];
        y_k=[y(i,j),y(i+1,j),y(i+1,j+1),y(i,j+1),y(i,j)];
        u_k=[u(i,j),u(i+1,j),u(i+1,j+1),u(i,j+1),u(i,j)];
        v_k=[v(i,j),v(i+1,j),v(i+1,j+1),v(i,j+1),v(i,j)];
        tmp_area=0;
        tmp_dudx=0;
        tmp_dudy=0;
        tmp_dvdx=0;
        tmp_dvdy=0;
        for k=1:4
            tmp_area=tmp_area+0.5*(x_k(k)*y_k(k+1)-y_k(k)*x_k(k+1));
            tmp_dudx=tmp_dudx+0.5*(u_k(k+1)+u_k(k))*(y_k(k+1)-y_k(k));
            tmp_dudy=tmp_dudy+0.5*(u_k(k+1)+u_k(k))*(x_k(k+1)-x_k(k));
            tmp_dvdx=tmp_dvdx+0.5*(v_k(k+1)+v_k(k))*(y_k(k+1)-y_k(k));
            tmp_dvdy=tmp_dvdy+0.5*(v_k(k+1)+v_k(k))*(x_k(k+1)-x_k(k));
        end
        dudx(i,j)=tmp_dudx/tmp_area;
        dudy(i,j)=-tmp_dudy/tmp_area;
        dvdx(i,j)=tmp_dvdx/tmp_area;
        dvdy(i,j)=-tmp_dvdy/tmp_area;
        tmp_area=abs(tmp_area);
        area(i,j)=tmp_area;
        scale(i,j)=sqrt(tmp_area);
        tri_regular_1(i,j)=(i  )+nb_x*(j  -1);
        tri_regular_2(i,j)=(i+1)+nb_x*(j  -1);
        tri_regular_3(i,j)=(i+1)+nb_x*(j+1-1);
        tri_regular_4(i,j)=(i  )+nb_x*(j+1-1);
    end
end

dudx=dudx(:);
dudy=dudy(:);
dvdx=dvdx(:);
dvdy=dvdy(:);
scale=scale(:);
area=area(:);

tri_regular(:,1)=tri_regular_1(:);
tri_regular(:,2)=tri_regular_2(:);
tri_regular(:,3)=tri_regular_3(:);
tri_regular(:,4)=tri_regular_4(:);

end
