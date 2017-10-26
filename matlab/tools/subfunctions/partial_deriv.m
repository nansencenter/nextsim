function [dudx,dudy,dvdx,dvdy,scale,area,minang,long_side]=partial_deriv(xy,uv)
%  function [dudx,dudy,dvdx,dvdy,scale,area]=partial_deriv(xy,uv)
%  Authors: Sylvain Bouillon, summer 2013
%           Lucas Girard, spring 2008
%  
%  GOAL: to compute partial derivatives, using the contour integral method (RGPS user's manual, R. Kwok)
%  
%  INPUT:
%  xy: coordinates of the 3 corners of each triangle
%  uv: velocity at the 3 corners of each triangle
%  
%  OUTPUT:
%  dudx, dudy,...: partial derivatives
%  scale, area : scale (=sqrt(area)) and area of the elements 

a=[xy(:,1,1),xy(:,1,2),0*xy(:,1,1)];
b=[xy(:,2,1),xy(:,2,2),0*xy(:,1,1)];
c=[xy(:,3,1),xy(:,3,2),0*xy(:,1,1)];
ab=[xy(:,2,:)-xy(:,1,:)];
ac=[xy(:,3,:)-xy(:,1,:)];
% cross product
ar=ab(:,1).*ac(:,2)-ac(:,1).*ab(:,2);
%ar=sum(cross(b-a,c-a,2),2);

lx=[ xy(:,1,1), xy(:,2,1), xy(:,3,1)];
ly=[ xy(:,1,2), xy(:,2,2), xy(:,3,2)];
lu=[ uv(:,1,1), uv(:,2,1), uv(:,3,1)];
lv=[ uv(:,1,2), uv(:,2,2), uv(:,3,2)];

% from Gunnar
% sides of triangle
ta = sqrt((xy(:,2,1)-xy(:,1,1)).^2+(xy(:,2,2)-xy(:,1,2)).^2);
tb = sqrt((xy(:,3,1)-xy(:,2,1)).^2+(xy(:,3,2)-xy(:,2,2)).^2);
tc = sqrt((xy(:,3,1)-xy(:,1,1)).^2+(xy(:,3,2)-xy(:,1,2)).^2);
% smallest angle
sides = [ta,tb,tc];
long_side = max(sides,[],2);
[shortest, id] = min(sides,[],2);% shortest side
for i=1:length(id)
    for k=1:3,
        if k~=id(i)
           sides(i,id(i)) = sides(i,k); % two other sides
        end
    end
end
sides=sides(:,1:2);
minang = acosd((sides(:,1).^2+sides(:,2).^2-shortest.^2)./2./sides(:,1)./sides(:,2)); % law of cosines
% end from Gunnar

% to order the coordinates in the counter-clockwise sense
ind=find(ar<0);
if(length(ind)>0)
    %ACB = counter-clockwise
    lx(ind,:)=[ xy(ind,1,1), xy(ind,3,1), xy(ind,2,1)];
    ly(ind,:)=[ xy(ind,1,2), xy(ind,3,2), xy(ind,2,2)];
    lu(ind,:)=[ uv(ind,1,1), uv(ind,3,1), uv(ind,2,1)];
    lv(ind,:)=[ uv(ind,1,2), uv(ind,3,2), uv(ind,2,2)];
end

% compute partial derivatives, area and scale for all the elements at once
area(:,1)=abs(ar)/2;
scale = sqrt(area);
[dudx(:,1),dudy(:,1),dvdx(:,1),dvdy(:,1)]=sub_diff_contour(lx,ly,lu,lv,abs(ar)/2);


function [dudx,dudy,dvdx,dvdy]=sub_diff_contour(x,y,u,v,aire)
%function [dudx,dudy,dvdx,dvdy]=sub_diff_contour_x,y,u,v,aire)
%compute partial derivates (x et y coordinates have to be ordered in the counter-clockwise sense)

temp = 1./(2.*aire);

dudx=  temp.*( (u(:,2)+u(:,1)).*(y(:,2)-y(:,1)) + (u(:,3)+u(:,2)).*(y(:,3)-y(:,2)) + (u(:,1)+u(:,3)).*(y(:,1)-y(:,3)) );
dudy= -temp.*( (u(:,2)+u(:,1)).*(x(:,2)-x(:,1)) + (u(:,3)+u(:,2)).*(x(:,3)-x(:,2)) + (u(:,1)+u(:,3)).*(x(:,1)-x(:,3)) );

dvdx=  temp.*( (v(:,2)+v(:,1)).*(y(:,2)-y(:,1)) + (v(:,3)+v(:,2)).*(y(:,3)-y(:,2)) + (v(:,1)+v(:,3)).*(y(:,1)-y(:,3)) );
dvdy= -temp.*( (v(:,2)+v(:,1)).*(x(:,2)-x(:,1)) + (v(:,3)+v(:,2)).*(x(:,3)-x(:,2)) + (v(:,1)+v(:,3)).*(x(:,1)-x(:,3)) );

return