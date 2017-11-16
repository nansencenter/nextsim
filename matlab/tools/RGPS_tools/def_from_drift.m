function [ntri,ttdnum,ttdeltat,area,xy_tricorner,uv,uv_tricorner,dudx,dudy,dvdx,dvdy,id_stream,quality_index]=def_from_drift(txy2,tdeltat2,tdnum2,trux2,trvy2,stream,min_area,max_area,min_minang,max_long_side,showplot,filter_noise,min_def,max_level)

% Triangulation
dt = DelaunayTri(txy2);
if showplot
    disp('Push button to continue!')
    pause
    close all
    figure; clf
    triplot(dt); axis image
end

ntri = length(dt(:,1));

% initialization
ttdnum = tdnum2;
ttdeltat = tdeltat2;
area=[];
dudx=[];
dudy=[];
dvdx=[];
dvdy=[];
xy_tricorner = zeros(ntri,3,2);
uv_tricorner = zeros(ntri,3,2);
uv = zeros(ntri,2);
id_stream = char(ones(ntri,1)*'z');
quality_index=0;

% Formatting the data
x=dt.X(:,1);
y=dt.X(:,2);
xy_tricorner(:,:,1) = x(dt(:,:)); % in m
xy_tricorner(:,:,2) = y(dt(:,:)); % in m    

uv_tricorner(:,:,1) = trux2(dt(:,:)); %m/day
uv_tricorner(:,:,2) = trvy2(dt(:,:)); %m/day
uv(:,:) = mean(uv_tricorner(:,:,:),2); %m/day

id_stream(:) = stream;

% compute derivatives and some properties of the triangles
[dudx,dudy,dvdx,dvdy,scale,area,minang,long_side]=partial_deriv(xy_tricorner,uv_tricorner);

%---------------- first selection ------------------------

% remove too big, too small triangles, and acute triangles
goodtri = find(area > min_area & area < max_area & minang > min_minang & long_side < max_long_side);
if(isempty(goodtri))
    ntri=-1;
    return
end

warning('OFF','MATLAB:TriRep:PtsNotInTriWarnId')
dt2 = TriRep(dt.Triangulation(goodtri,:),dt.X);
warning('ON','MATLAB:TriRep:PtsNotInTriWarnId')

% reduced mesh
tri=dt2.Triangulation;
x=dt2.X(:,1);
y=dt2.X(:,2);
warning('OFF')
tr=TriRep(tri,x(:),y(:));
warning('ON')

%---------------- second selection ------------------------

nb_triangles=size(tri,1);

% We search the triangles or pairs of triangles that are alone 
triangles_alone=zeros(nb_triangles,1);
next_triangle=neighbors(tr,(1:nb_triangles)');
number_no_neighbors=sum(isnan(next_triangle),2);

% a triangle alone in the blue (not taken into account)
triangles_alone(number_no_neighbors==3,1)=1;

% a couple of triangles alone in the blue (not taken into account)
ind_candidate=find(number_no_neighbors==2);
next_triangle_tmp=next_triangle(ind_candidate,:);
next_triangle_tmp=next_triangle_tmp(~isnan(next_triangle_tmp));
[indice_double] = intersect(ind_candidate,next_triangle_tmp);
triangles_alone(indice_double,1)=1;

% only keep the not alone elements
not_alone=find(triangles_alone==0);
if(isempty(not_alone))
    ntri=-2;
    return
end

warning('OFF','MATLAB:TriRep:PtsNotInTriWarnId')
dt3 = TriRep(dt2.Triangulation(not_alone,:),dt2.X);
warning('ON','MATLAB:TriRep:PtsNotInTriWarnId')
dt2=dt3;

tri=dt2.Triangulation;
x=dt2.X(:,1);
y=dt2.X(:,2);
warning('OFF')
tr=TriRep(tri,x(:),y(:));
warning('ON')

%---------------- save the selected triangles ------------------------

ntri = length(goodtri(not_alone));

area = area(goodtri(not_alone));
minang = minang(goodtri(not_alone));
dudx = dudx(goodtri(not_alone)); dudy = dudy(goodtri(not_alone));
dvdx = dvdx(goodtri(not_alone)); dvdy = dvdy(goodtri(not_alone));
xy_tricorner = xy_tricorner(goodtri(not_alone),:,:);
uv_tricorner = uv_tricorner(goodtri(not_alone),:,:);
uv = uv(goodtri(not_alone),:);
id_stream = id_stream(goodtri(not_alone));

% compute the invariants and the deformation criterion for
% smoothing
[invar]=invariants(dudx(:,1),dudy(:,1),dvdx(:,1),dvdy(:,1));

    
def=invar.eps';

deformed=find(def>=min_def);
is_deformed=zeros(size(invar.div));
is_deformed(deformed)=1;

if(showplot && (ntri > 3))
    RGPS_plots_veloc(dt2,uv);
    RGPS_plots_deriv(dt2,invar,deformed,1)
end % /showplot

if(filter_noise)
    
    %-------------------SMART SMOOTHING---------------------%
 
    new_dudx=dudx;
    new_dudy=dudy;
    new_dvdx=dvdx;
    new_dvdy=dvdy;
    
    nb_other_deformed=zeros(size(invar.div));
    
    for i=1:length(deformed)
        
        other_deformed=recursive_detect_feature([],deformed(i),1,dt2,is_deformed,max_level);
        other_deformed=unique(other_deformed);
        nb_other_deformed(deformed(i))=length(other_deformed);
        
        if(~isempty(other_deformed))
            
            area_other_deformed=sum(area(other_deformed));
            
            new_dudx(deformed(i),1)=(dudx(other_deformed)'*area(other_deformed))/area_other_deformed;
            new_dudy(deformed(i),1)=(dudy(other_deformed)'*area(other_deformed))/area_other_deformed;
            new_dvdx(deformed(i),1)=(dvdx(other_deformed)'*area(other_deformed))/area_other_deformed;
            new_dvdy(deformed(i),1)=(dvdy(other_deformed)'*area(other_deformed))/area_other_deformed;
        end
    end
    
    quality_index=sum((nb_other_deformed>=(max_level+1)).*(nb_other_deformed<=(4*max_level+1)))/sum(nb_other_deformed>0);
    
    [invar]=invariants(new_dudx,new_dudy,new_dvdx,new_dvdy);
    
    if(showplot && (ntri > 3))
        disp(quality_index)
        RGPS_plots_deriv(dt2,invar,deformed,2)
    end % /showplot
    
    dudx=new_dudx;
    dudy=new_dudy;
    dvdx=new_dvdx;
    dvdy=new_dvdy;
    
else
    quality_index=0;
end % /filter_noise


end

function [result]=recursive_detect_feature(visited_triangle,triangle,level,tr,is_deformed,max_level)

result=triangle;
if(level>max_level)
    return;
else
    next_triangle=neighbors(tr,triangle);
    % % check (it is commented as it takes a bit of time)
    % if(sum(isnan(next_triangle))==3) % a triangle alone in the blue (not taken into account)
    %     disp(tr)
    %     disp(next_triangle)
    %     disp(triangle)
    %     error('still a triangle alone in the blue !!!');
    % end
    
    visited_triangle=[visited_triangle,triangle];
    for k=1:3,
        next_ti=next_triangle(k);
        if(~isnan(next_ti))
            if(sum(next_ti==visited_triangle)==0)
                if(is_deformed(next_ti))
                    [tmp_result]=recursive_detect_feature(visited_triangle,next_ti,level+1,tr,is_deformed,max_level);
                    result=[result,tmp_result];
                end
            end
        end
    end
end
end

function RGPS_plots_veloc(dt2,uv)

allx1 = dt2.X(dt2.Triangulation(:,1),1);
allx2 = dt2.X(dt2.Triangulation(:,2),1);
allx3 = dt2.X(dt2.Triangulation(:,3),1);
allx = [allx1, allx2, allx3];
ally1 = dt2.X(dt2.Triangulation(:,1),2);
ally2 = dt2.X(dt2.Triangulation(:,2),2);
ally3 = dt2.X(dt2.Triangulation(:,3),2);
ally = [ally1, ally2, ally3];

% x-component of the velocity anomaly

figure(15); clf
fill(allx',ally',(uv(:,1)'-mean(uv(:,1)))/1000,'EdgeColor','none')
caxis([-1.5 1.5]);
axis image
colormap('blue2red');
t=colorbar;
set(get(t,'ylabel'),'String', '[km/day]');
title('x-component of the velocity anomaly')
xlabel('x [km]')
ylabel('y [km]')

axis equal

% y-component of the velocity anomaly

figure(16); clf
fill(allx',ally',(uv(:,2)'-mean(uv(:,2)))/1000,'EdgeColor','none')
caxis([-1.5 1.5]);
axis image
colormap('blue2red');
t=colorbar;
set(get(t,'ylabel'),'String', '[km/day]');
title('y-component of the velocity anomaly')
xlabel('x [km]')
ylabel('y [km]')
axis equal
end

function RGPS_plots_deriv(dt2,invar,deformed,plot_deformed)
allx1 = dt2.X(dt2.Triangulation(:,1),1);
allx2 = dt2.X(dt2.Triangulation(:,2),1);
allx3 = dt2.X(dt2.Triangulation(:,3),1);
allx = [allx1, allx2, allx3];
ally1 = dt2.X(dt2.Triangulation(:,1),2);
ally2 = dt2.X(dt2.Triangulation(:,2),2);
ally3 = dt2.X(dt2.Triangulation(:,3),2);
ally = [ally1, ally2, ally3];

% figure(32); clf
% triplot(dt2); axis image
% [fe xe] = freeBoundary(dt);
% hold on
% plot(xe(:,1), xe(:,2), '-r', 'LineWidth',2) ;

% divergence plot
figure(plot_deformed*10+1);
fill(allx',ally',invar.div','EdgeColor','none')
if(plot_deformed==3)
    hold on
end
caxis([-0.04 0.04]); axis image
colormap('blue2red');
t=colorbar;
set(get(t,'ylabel'),'String', '[/day]');
title('Divergence rate')
xlabel('x [km]')
ylabel('y [km]')
axis equal

% shear plot
figure(plot_deformed*10+2);
fill(allx',ally',invar.shear','EdgeColor','none')
if(plot_deformed==3)
    hold on
end
caxis([0 0.08]); axis image
colormap('gray2red');
t=colorbar;
set(get(t,'ylabel'),'String', '[/day]');
title('Shear rate')
xlabel('x [km]')
ylabel('y [km]')
axis equal

% shear plot
figure(plot_deformed*10+3);
fill(allx',ally',invar.eps','EdgeColor','none')
if(plot_deformed==1)
    hold on
    fill(allx(deformed,:)',ally(deformed,:)',invar.eps(deformed,:)','EdgeColor','black')
    hold off
end
if(plot_deformed==3)
    hold on
end
caxis([0 0.08]); axis image
colormap('gray2red');
t=colorbar;
set(get(t,'ylabel'),'String', '[/day]');
title('Total deformation rate')
xlabel('x [km]')
ylabel('y [km]')
axis equal

%     % vorticity plot
%
%     figure;
%     fill(allx',ally',invar.vor','EdgeColor','none')
%     caxis([-0.04 0.04]);
%     axis image
%     colormap('blue2red');
%     t=colorbar;
%     set(get(t,'ylabel'),'String', '[/day]');
%     title('vorticity')
%     xlabel('x [km]')
%     ylabel('y [km]')
%     axis equal

end
