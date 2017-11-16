
function test_div()
%close all
clear all
%figure

% ---- Parameters
showplot=1;
resolution=0.1;
max_nb_test=0;
min_def=1e-9;
nb_of_fractures=1;
scaling_displ=1;
test_ok=0;
% ---- end Parameters


if(showplot)
    nb_test_i=1;
    vec_max_level=[0:25];%0:7]
else
    nb_test_i=100;
    vec_max_level=[25];
end

error_div_no_swap=zeros(nb_test_i,length(vec_max_level));
error_conv_no_swap=zeros(nb_test_i,length(vec_max_level));
error_div=zeros(nb_test_i,length(vec_max_level));

k=1;
for test_i=1:nb_test_i,
    
    if(nb_of_fractures==1)
        dx=[0.01];
        dy=[0.00];
    elseif(nb_of_fractures==2)
        dx_tmp=0.01;
        dy_tmp=0; %dx_tmp/8;
        
        dx=[dx_tmp,-dx_tmp];
        
        %dx=[dx_tmp,dy_tmp];
        dy=[-dy_tmp,0.00];
    else
        error('nb of fractures ~=1 or 2 is not coded')
    end

    [x,y] = meshgrid(0:resolution:1,0:resolution:1);
    [dim1,dim2]=size(x);
    x(2:dim1-1,2:dim2-1)=x(2:dim1-1,2:dim2-1)+(rand(size(x(2:dim1-1,2:dim2-1)))-0.5)*2*resolution/3;
    y(2:dim1-1,2:dim2-1)=y(2:dim1-1,2:dim2-1)+(rand(size(y(2:dim1-1,2:dim2-1)))-0.5)*2*resolution/3;
    
    u=0*x;
    v=0*x;

    disp(test_i);
    
    xc(1)=0.5;
    yc(1)=0.5;
    if(    test_ok==0)
        angle_i(1)=0;
    elseif(test_ok==1)
        angle_i(1)=atan(0.15/1);
    elseif(test_ok==3)
        angle_i(1)=atan(0.2/1)-1e-8;
    else
        angle_i(1)=(-1+2*rand)*(atan(0.2/1)-1e-8);
    end
    
    
    if(length(dx)>1)
        xc(2)=0.5;
        yc(2)=0.5;
        angle_i(2)=atan(0.1/1)+rand*(pi/2);
    end
    
    for i=1:length(dx),

        Dx=resolution;
        theta=angle_i(i);
        L=1;
        l =sqrt(L^2+(L*tan(theta)   )^2);
        lp=sqrt(L^2+(L*tan(theta)-Dx)^2);
        dtheta=theta-atan(tan(theta)-Dx/L);
        error_div=lp*dx*sin(dtheta);
        approm_error=dx*Dx;
        
        alpha(i)=tan(angle_i(i));
        
        d(1)=(1-xc(i))/cos(angle_i(i));
        d(2)=(1-yc(i))/sin(angle_i(i));
        d(3)=(0-xc(i))/cos(angle_i(i));
        d(4)=(0-yc(i))/sin(angle_i(i));
        
        length_fract(i)=min(d(d>0))-max(d(d<0));
        
        %         xc = 0.5784;
        %         yc = 0.6311;
        %         alpha = -1.6766;
        
        % test bug concavity
        % xc =    [0.6309    0.6437];
        % yc =    [0.5815    0.6937];
        % alpha =    [1.5755    0.0988];
        
        if(i>1)
            if(angle_i(2)<pi/2)
                ind=find((((y-yc(i))-alpha(i)*(x-xc(i)))<0).*((y-yc(1))-alpha(1)*(x-xc(1)))>0);
            else
                ind=find((((y-yc(i))-alpha(i)*(x-xc(i)))>0).*((y-yc(1))-alpha(1)*(x-xc(1)))>0);
            end
        else
            ind=find(((y-yc(i))-alpha(i)*(x-xc(i)))>0);
        end
        
%         u(ind)=u(ind)+(dx(i)*cos(atan(alpha(i)))-dy(i)*sin(atan(alpha(i))));
%         v(ind)=v(ind)+(dx(i)*sin(atan(alpha(i)))+dy(i)*cos(atan(alpha(i))));
        u(ind)=u(ind)+(dx(i)*cos(atan(alpha(1)))-dy(i)*sin(atan(alpha(1))));
        v(ind)=v(ind)+(dx(i)*sin(atan(alpha(1)))+dy(i)*cos(atan(alpha(1))));
    end
    
    tri=delaunay(x(:),y(:));
    tr=TriRep(tri,x(:),y(:));
    e=edges(tr);
    ti = edgeAttachments(tr,e);
    
    % deformation with the triangulation
    xy_tricorner(:,:,1)=x(tri);
    xy_tricorner(:,:,2)=y(tri);
    uv_tricorner(:,:,1)=u(tri);
    uv_tricorner(:,:,2)=v(tri);
    
    x_center=mean(xy_tricorner(:,:,1),2);
    y_center=mean(xy_tricorner(:,:,2),2);
    
    [data.dudx,data.dudy,data.dvdx,data.dvdy,data.scale,data.area]=partial_deriv(xy_tricorner,uv_tricorner);
    [invar]=invariants(data.dudx,data.dudy,data.dvdx,data.dvdy);
       
    if(showplot)
        [k]=plots(tri,invar,angle_i,xc,yc,x,u,y,v,dx,dy,k,scaling_displ,[],[]);
    end
       
    [error_div_no_swap_i,error_conv_no_swap_i,k,ex_deformed]=smart_smoothing(showplot,vec_max_level,data,invar,tr,x_center,y_center,tri,angle_i,xc,yc,x,u,y,v,dx,dy,k,scaling_displ,length_fract,min_def,-1);
    error_div_no_swap(test_i,:)=error_div_no_swap_i;
    error_conv_no_swap(test_i,:)=error_conv_no_swap_i;
   
    [data_regular.dudx,data_regular.dudy,data_regular.dvdx,data_regular.dvdy,data_regular.scale,data_regular.area,tri_regular]=partial_deriv_regular(x,y,u,v);
    [invar_regular]=invariants(data_regular.dudx,data_regular.dudy,data_regular.dvdx,data_regular.dvdy);
    
    if(length(dx)==1)
        error_div_regular (test_i)=sum(invar_regular.div(find(invar_regular.div>0)).*data_regular.area(find(invar_regular.div>0))) -length_fract(1)*dy(1)*(dy(1)>0);
        error_conv_regular(test_i)=sum(invar_regular.div(find(invar_regular.div<0)).*data_regular.area(find(invar_regular.div<0)))-length_fract(1)*dy(1)*(dy(1)<0);   
        error_div_regular (test_i)=sum(invar_regular.div(find(invar_regular.div>0)).*data_regular.area(find(invar_regular.div>0))) -length_fract(1)*dy(1)*(dy(1)>0);
        error_conv_regular(test_i)=sum(invar_regular.div(find(invar_regular.div<0)).*data_regular.area(find(invar_regular.div<0)))-length_fract(1)*dy(1)*(dy(1)<0);   
    else
        error_div_regular (test_i)=sum(invar_regular.div(find(invar_regular.div>0)).*data_regular.area(find(invar_regular.div>0))) -(length_fract(1)*dy(1)*(dy(1)>0)+length_fract(2)/2*dx(2)*(dx(2)>0)*sin(angle_i(2)-angle_i(1)));
        error_conv_regular(test_i)=sum(invar_regular.div(find(invar_regular.div<0)).*data_regular.area(find(invar_regular.div<0)))-(length_fract(1)*dy(1)*(dy(1)<0)+length_fract(2)/2*dx(2)*(dx(2)<0)*sin(angle_i(2)-angle_i(1)));
    end
        
%     if(showplot)
%         [error_div_no_swap_no_detect_i,error_conv_no_swap_no_detect_i,k,ex_deformed]=smart_smoothing(showplot,vec_max_level,data,invar,tr,x_center,y_center,tri,angle_i,xc,yc,x,u,y,v,dx,dy,k,scaling_displ,length_fract,0,ex_deformed);
%         error_div_no_swap_no_detect(test_i,:)=error_div_no_swap_i;
%         error_conv_no_swap_no_detect(test_i,:)=error_conv_no_swap_i;
%         
%         plots(tri_regular,invar_regular,angle_i,xc,yc,x,u,y,v,dx,dy,k,scaling_displ,[],[]);
%     end

    if(max_nb_test>0)
        
        nb_test=0;
        test_again=1;
        while(test_again && nb_test<max_nb_test)
            nb_test=nb_test+1;
            new_tri=tri;
            old_cells=[];
            ind=find(cellfun(@length, ti)==2);
            nb_ugly=0;
            
            cells_ind=cell2mat(ti(ind));
            mean_abs_div = mean(abs(invar.div(cells_ind)),2);
            [Y,i_sort] = sort(mean_abs_div,'descend');
            
            for i=1:length(ind),
                i_test=i_sort(i);
                %cells_id=cell2mat(ti(ind(i)));
                cells_id=cells_ind(i_test,:);
                if(~length(find(cells_id(1)==old_cells)))
                    if(~length(find(cells_id(2)==old_cells)))
                        
                        
                        
                        max_div=max(abs(invar.div(cells_id)));
                        mean_div=mean(abs(invar.div(cells_id)));
                        %if((invar.div(cells_id(1))*invar.div(cells_id(2)))<0)
                        if(max_div>1e-7)
                            
                            ind1=find((tri(cells_id(1),:)~=e(ind(i_test),1)).*(tri(cells_id(1),:)~=e(ind(i_test),2)));
                            ind2=find((tri(cells_id(2),:)~=e(ind(i_test),1)).*(tri(cells_id(2),:)~=e(ind(i_test),2)));
                            
                            e1=e(ind(i_test),1); % first node of the common edge
                            e2=e(ind(i_test),2); % second node of the common edge
                            c1=tri(cells_id(1),ind1); % third node for cell 1
                            c2=tri(cells_id(2),ind2); % third node for cell 2
                            
                            tri_tmp(1,:)=[c1,c2,e1];
                            tri_tmp(2,:)=[c1,c2,e2];
                            
                            xy_tricorner_tmp(:,:,1)=x(tri_tmp);
                            xy_tricorner_tmp(:,:,2)=y(tri_tmp);
                            uv_tricorner_tmp(:,:,1)=u(tri_tmp);
                            uv_tricorner_tmp(:,:,2)=v(tri_tmp);
                            
                            [dudx,dudy,dvdx,dvdy,scale,area]=partial_deriv(xy_tricorner_tmp,uv_tricorner_tmp);
                            [invar_tmp]=invariants(dudx,dudy,dvdx,dvdy);
                            
                            % we must check if the two previous triangles formed a convex quadrangle (otherwise the swap is ill defined)
                            % Compute the inner angle at the point e1
                            vec_e1e2=[x(e2)-x(e1),y(e2)-y(e1)];
                            vec_e1c1=[x(c1)-x(e1),y(c1)-y(e1)];
                            vec_e1c2=[x(c2)-x(e1),y(c2)-y(e1)];
                            
                            angle_e1c1_on_e1e2=acos(dot(vec_e1c1,vec_e1e2)/(norm(vec_e1c1)*norm(vec_e1e2)));
                            angle_e1c2_on_e1e2=acos(dot(vec_e1c2,vec_e1e2)/(norm(vec_e1c2)*norm(vec_e1e2)));
                            angle_interior_1 = abs(angle_e1c1_on_e1e2)+abs(angle_e1c2_on_e1e2);
                            
                            % Compute the inner angle at the point e2
                            vec_e2e1=[x(e1)-x(e2),y(e1)-y(e2)];
                            vec_e2c1=[x(c1)-x(e2),y(c1)-y(e2)];
                            vec_e2c2=[x(c2)-x(e2),y(c2)-y(e2)];
                            
                            angle_e2c1_on_e2e1=acos(dot(vec_e2c1,vec_e2e1)/(norm(vec_e2c1)*norm(vec_e2e1)));
                            angle_e2c2_on_e2e1=acos(dot(vec_e2c2,vec_e2e1)/(norm(vec_e2c2)*norm(vec_e2e1)));
                            angle_interior_2 = abs(angle_e2c1_on_e2e1)+abs(angle_e2c2_on_e2e1);
                            
                            %if ( min(minang)>10 && max(abs(invar_tmp.div))<max_div )
                            if ( min(minang)>10 && mean(abs(invar_tmp.div))<0.9*mean_div && angle_interior_1<pi && angle_interior_2<pi)
                                nb_ugly=nb_ugly+1;
                                new_tri(cells_id(1),:)=tri_tmp(1,:);
                                new_tri(cells_id(2),:)=tri_tmp(2,:);
                                old_cells=[old_cells,cells_id(1),cells_id(2)];
                            end
                        end
                    end
                end
            end
            nb_ugly
            if(nb_ugly==0)
                test_again=0;
            else
                
                tri=new_tri;
                
                tr=TriRep(tri,x(:),y(:));
                e=edges(tr);
                ti = edgeAttachments(tr,e);
                
                xy_tricorner(:,:,1)=x(tri);
                xy_tricorner(:,:,2)=y(tri);
                uv_tricorner(:,:,1)=u(tri);
                uv_tricorner(:,:,2)=v(tri);
                
                [data.dudx,data.dudy,data.dvdx,data.dvdy,data.scale,data.area]=partial_deriv(xy_tricorner,uv_tricorner);
                [invar]=invariants(data.dudx,data.dudy,data.dvdx,data.dvdy);
                
            end
        end
                
        plots(tri,invar,angle_i,xc,yc,x,u,y,v,dx,dy,k,scaling_displ,[],[])
        
        [error_div,error_conv,k,ex_defomed]=smart_smoothing(showplot,vec_max_level,data,invar,tr,x_center,y_center,tri,angle_i,xc,yc,x,u,y,v,dx,dy,k,scaling_displ,length_fract,min_def,-1);
        error_div(test_i,:)=error_div_i;
        error_conv(test_i,:)=error_conv_i;
    end
    
end


for j=1:length(vec_max_level)
    if(max_nb_test>0)
        median_error_div(j)=median(error_div(:,j));
        median_error_conv(j)=median(error_conv(:,j));
        rms_error_div(j)=rms(error_div(:,j));
        rms_error_conv(j)=rms(error_conv(:,j));
    end

    median_error_div_no_swap(j)=median(error_div_no_swap(:,j));
    median_error_conv_no_swap(j)=median(error_conv_no_swap(:,j));
    rms_error_div_no_swap(j)=rms(error_div_no_swap(:,j));
    rms_error_conv_no_swap(j)=rms(error_conv_no_swap(:,j));
end

median_error_div_regular=median(error_div_regular);
median_error_conv_regular=median(error_conv_regular);

rms_error_div_regular=rms(error_div_regular);
rms_error_conv_regular=rms(error_conv_regular);

% figure
% hold on
% for test_i=1:nb_test_i,
%     if(max_nb_test>0)        
%         plot(vec_max_level,error_div(test_i,:))
%     end
%     plot(vec_max_level,error_div_no_swap(test_i,:))
%     
% end
% if(max_nb_test>0)
%     
%     plot(vec_max_level,median_error_div,'r','LineWidth',4)
% end
% plot(vec_max_level,median_error_div_no_swap,'g','LineWidth',4)

% figure
% hold on
% for test_i=1:nb_test_i,
%     if(max_nb_test>0)        
%         plot(vec_max_level,error_conv(test_i,:))
%     end
%     plot(vec_max_level,error_conv_no_swap(test_i,:))
%     
% end



figure
hold on

error_tot_no_swap=abs(error_conv_no_swap)+abs(error_div_no_swap);
error_tot_regular=abs(error_conv_regular)+abs(error_div_regular);
bias_tot_no_swap=error_conv_no_swap+error_div_no_swap;
bias_tot_regular=error_conv_regular+error_div_regular;

% for test_i=1:nb_test_i,
%     plot(vec_max_level,error_tot_no_swap(test_i,:));
%     plot(0,error_tot_regular(test_i),'sr','MarkerFaceColor','r','MarkerSize',4);
% end

median_error_no_swap=median(error_tot_no_swap,1);
median_error_regular=median(error_tot_regular);

rms_error_no_swap=rms(error_tot_no_swap,1);
rms_error_regular=rms(error_tot_regular);
plot(vec_max_level,median_error_no_swap,'g','LineWidth',4)
plot(0,median_error_regular,'sr','MarkerFaceColor','r','MarkerSize',10)

x_axis=xlim;
y_axis=ylim;
axis([x_axis 0 y_axis(2)])

figure
hold on

plot(vec_max_level,rms_error_no_swap,'g','LineWidth',4)
plot(0,rms_error_regular,'sr','MarkerFaceColor','r','MarkerSize',10)

x_axis=xlim;
y_axis=ylim;
axis([x_axis 0 y_axis(2)])

% if(max_nb_test>0)
%     
%     plot(vec_max_level,median_error_conv,'r','LineWidth',4)
% end
% plot(vec_max_level,median_error_conv_no_swap,'g','LineWidth',4)

if(showplot)
    
    
    disp('bias_tot_regular');
    bias_tot_regular
    disp('error_tot_regular');
    error_tot_regular
    disp('median_error_regular');
    median_error_regular
    disp('error_div_regular');
    error_div_regular
    disp('median_error_div_regular');
    median_error_div_regular
    disp('error_conv_regular');
    error_conv_regular
    disp('median_error_conv_regular')
    median_error_conv_regular
    
    disp('bias_tot_no_swap');
    bias_tot_no_swap
    disp('error_tot_no_swap');
    error_tot_no_swap
    disp('median_error_no_swap');
    median_error_no_swap
    disp('vec_max_level');
    error_div_no_swap
    disp('median_error_div_no_swap');
    median_error_div_no_swap
    disp('error_conv_no_swap');
    error_conv_no_swap
    disp('median_error_conv_no_swap')
    median_error_conv_no_swap
    
    vec_max_level
    disp('error_div_no_swap');
    
    else
    save(['error_',num2str(nb_of_fractures),'_sliding_frac_dx_01_dy_00_',num2str(test_ok),'_',num2str(1/resolution)],'bias_tot_no_swap','bias_tot_regular','error_tot_regular','median_error_regular','error_tot_no_swap','median_error_no_swap','vec_max_level','error_div_no_swap','median_error_div_no_swap','error_conv_no_swap','median_error_conv_no_swap')
end

error('end')


%------------ third test correction (split bad elements to minimize the absolute divergence) --------------

data_tris=data;

nb_init=length(tri(:,1));
nb_deformed=length(deformed);

% we will split all the deformed triangle in two so we need nb_deformed new
% element
tri_new=[tri;zeros(nb_deformed,3)];
x_new=[x(:);zeros(nb_deformed,1)];
u_new=[u(:);zeros(nb_deformed,1)];
y_new=[y(:);zeros(nb_deformed,1)];
v_new=[v(:);zeros(nb_deformed,1)];

for i=1:nb_deformed,
    
    div_min_max=1e9;%invar.div(deformed(i));
    
    j_other_best=0;
    j_best=0;
    
    ind_node=[1:3];
    ind_next_node=mod([1:3]-1+1,3)+1;
    ind_prev_node=mod([1:3]-1-1,3)+1;
    
    dx_to_next=xy_tricorner(deformed(i),ind_next_node,1)-xy_tricorner(deformed(i),ind_node,1);
    dy_to_next=xy_tricorner(deformed(i),ind_next_node,2)-xy_tricorner(deformed(i),ind_node,2);
    dx_to_prev=xy_tricorner(deformed(i),ind_prev_node,1)-xy_tricorner(deformed(i),ind_node,1);
    dy_to_prev=xy_tricorner(deformed(i),ind_prev_node,2)-xy_tricorner(deformed(i),ind_node,2);
    
    orientation_to_next=atan2(dy_to_next,dx_to_next);
    orientation_to_prev=atan2(dy_to_prev,dx_to_prev);
    
    diff_orient_to_prev=mod(orientation_to_prev-orientation(i),2*pi)
    diff_orient_to_next=mod(orientation(i)-orientation_to_next,2*pi)
    diff_orient_to=mod(orientation_to_prev-orientation_to_next,2*pi)
    
    dx_from_next=xy_tricorner(deformed(i),ind_node,1)-xy_tricorner(deformed(i),ind_next_node,1);
    dy_from_next=xy_tricorner(deformed(i),ind_node,2)-xy_tricorner(deformed(i),ind_next_node,2);
    dx_from_prev=xy_tricorner(deformed(i),ind_node,1)-xy_tricorner(deformed(i),ind_prev_node,1);
    dy_from_prev=xy_tricorner(deformed(i),ind_node,2)-xy_tricorner(deformed(i),ind_prev_node,2);
    
    orientation_from_next=atan2(dy_from_next,dx_from_next);
    orientation_from_prev=atan2(dy_from_prev,dx_from_prev);
    
    diff_orient_from_prev=mod(orientation_from_prev-orientation(i),2*pi)
    diff_orient_from_next=mod(orientation(i)-orientation_from_next,2*pi)
    diff_orient_from=mod(orientation_from_prev-orientation_from_next,2*pi)
    
    for j=1:4,
        if(j==4)
            error('bugbug')
        end
        if((diff_orient_from_prev(j)<=diff_orient_from(j)) && (diff_orient_from_next(j)<=diff_orient_from(j)))
            %             [min_diff,ind_min_diff]=min([diff_orient_from_prev(j),diff_orient_from_next(j)]);
            %             [max_diff,ind_max_diff]=max([diff_orient_from_prev(j),diff_orient_from_next(j)]);
            %
            %             k=[-1,+1];
            %             j_other=mod(j-1+k(ind_min_diff),3)+1
            %             j_other_end=mod(j-1+k(ind_max_diff),3)+1
            break
        end
        if((diff_orient_to_prev(j)  <=diff_orient_to(j)  ) && (diff_orient_to_next(j)  <=diff_orient_to(j)  ))
            %             [min_diff,ind_min_diff]=min([diff_orient_to_prev(j),diff_orient_to_next(j)]);
            %             [max_diff,ind_max_diff]=max([diff_orient_to_prev(j),diff_orient_to_next(j)]);
            %
            %             k=[-1,+1];
            %             j_other=mod(j-1+k(ind_min_diff),3)+1
            %             j_other_end=mod(j-1+k(ind_max_diff),3)+1
            break;
        end
    end
    
    
    % for each node
    %     for j=1:3,
    %
    %
    % %         d_theta_1   =mod(orientation(i)-orientation_edge(ind_node(j)),pi)
    %         d_theta_2   =mod(orientation(i)-orientation_edge(ind_prev_node(j)),pi)
    %         d_theta_1_2 =mod(orientation_edge(ind_node(j))-orientation_edge(ind_prev_node(j)),pi)
    %
    %         if(max(d_theta_1,d_theta_2)>d_theta_1_2 )
    %             if(j<3)
    %                 continue;
    %             else
    %                 error('error')
    %             end
    %         end
    %         i
    %         j
    %         orientation(i)
    
    %         d_theta_1   =mod(orientation(i)-orientation_edge(ind_node(j)),pi)
    %         d_theta_2   =mod(orientation(i)-orientation_edge(ind_prev_node(j)),pi)
    %
    %         k=[+1,-1];
    %         [min_k,ind_min_k]=min([d_theta_1,d_theta_2]);
    %         [max_k,ind_max_k]=max([d_theta_1,d_theta_2]);
    %
    %         j_other=mod(j-1+k(ind_min_k),3)+1
    %         j_other_end=mod(j-1+k(ind_max_k),3)+1
    
    %         dx_test=(1-beta)*x(tri(deformed(i),j_other))+beta*x(tri(deformed(i),j_other_end))-x(tri(deformed(i),j);
    %         dy_test=(1-beta)*y(tri(deformed(i),j_other))+beta*y(tri(deformed(i),j_other_end))-y(tri(deformed(i),j);
    %
    %         dx_test=beta*(x(tri(deformed(i),j_other_end))-x(tri(deformed(i),j_other)))+(x(tri(deformed(i),j_other))-x(tri(deformed(i),j));
    %         dy_test=beta*(y(tri(deformed(i),j_other_end))-y(tri(deformed(i),j_other)))+(y(tri(deformed(i),j_other))-y(tri(deformed(i),j));
    
    
    %         % we add a node on the edge towards the previous or next points
    %         for k=[-1,+1]
    %             j_other=mod(j-1+k,3)+1;
    %
    %             for beta=d_beta:d_beta:1-d_beta,
    %                 x_test=(1-beta)*x(tri(deformed(i),j))+beta*x(tri(deformed(i),j_other));
    %                 y_test=(1-beta)*y(tri(deformed(i),j))+beta*y(tri(deformed(i),j_other));
    %                 u_test=u(tri(deformed(i),j));
    %                 v_test=v(tri(deformed(i),j));
    %
    %                 xy_tricorner_bis(1,:,1)=x(tri(deformed(i),:));
    %                 xy_tricorner_bis(2,:,1)=x(tri(deformed(i),:));
    %                 xy_tricorner_bis(1,:,2)=y(tri(deformed(i),:));
    %                 xy_tricorner_bis(2,:,2)=y(tri(deformed(i),:));
    %
    %                 uv_tricorner_bis(1,:,1)=u(tri(deformed(i),:));
    %                 uv_tricorner_bis(2,:,1)=u(tri(deformed(i),:));
    %                 uv_tricorner_bis(1,:,2)=v(tri(deformed(i),:));
    %                 uv_tricorner_bis(2,:,2)=v(tri(deformed(i),:));
    %
    %                 xy_tricorner_bis(1,j_other,1)=x_test;
    %                 xy_tricorner_bis(2,j      ,1)=x_test;
    %                 xy_tricorner_bis(1,j_other,2)=y_test;
    %                 xy_tricorner_bis(2,j      ,2)=y_test;
    %
    %                 uv_tricorner_bis(1,j_other,1)=u_test;
    %                 uv_tricorner_bis(2,j      ,1)=u_test;
    %                 uv_tricorner_bis(1,j_other,2)=v_test;
    %                 uv_tricorner_bis(2,j      ,2)=v_test;
    %
    %                 [data_bis.dudx,data_bis.dudy,data_bis.dvdx,data_bis.dvdy,data_bis.scale,data_bis.area]=partial_deriv(xy_tricorner_bis,uv_tricorner_bis);
    %                 [invar_bis]=invariants(data_bis.dudx,data_bis.dudy,data_bis.dvdx,data_bis.dvdy);
    %
    %                 if(max(abs(invar_bis.div))<=abs(div_min_max))
    %                 %if(max(abs(invar_bis.eps))<abs(div_min_max))
    %                     div_min_max=max(abs(invar_bis.div));
    %                     %div_min_max=max(abs(invar_bis.eps));
    %
    %                     j_other_best=j_other;
    %                     j_best=j;
    %
    %                     x_best=x_test;
    %                     y_best=y_test;
    %                     u_best=u_test;
    %                     v_best=v_test;
    %                 else
    %                     break;
    %                 end
    %             end
    %        end
    
    for k=[-1,+1]
        j_other=mod(j-1+k,3)+1
        j_other_end=mod(j-1-k,3)+1
        
        a_y=(y(tri(deformed(i),j_other_end))-y(tri(deformed(i),j_other)));
        b_y=(y(tri(deformed(i),j_other))-y(tri(deformed(i),j)));
        
        a_x=(x(tri(deformed(i),j_other_end))-x(tri(deformed(i),j_other)));
        b_x=(x(tri(deformed(i),j_other))-x(tri(deformed(i),j)));
        
        c=tan(orientation(i));
        
        beta=-(b_y-c*b_x)/(a_y-c*a_x);
        %beta=0.5;
        
        x_test=(1-beta)*x(tri(deformed(i),j_other))+beta*x(tri(deformed(i),j_other_end));
        y_test=(1-beta)*y(tri(deformed(i),j_other))+beta*y(tri(deformed(i),j_other_end));
        u_test=u(tri(deformed(i),j_other));
        v_test=v(tri(deformed(i),j_other));
        
        xy_tricorner_bis(1,:,1)=x(tri(deformed(i),:));
        xy_tricorner_bis(2,:,1)=x(tri(deformed(i),:));
        xy_tricorner_bis(1,:,2)=y(tri(deformed(i),:));
        xy_tricorner_bis(2,:,2)=y(tri(deformed(i),:));
        
        uv_tricorner_bis(1,:,1)=u(tri(deformed(i),:));
        uv_tricorner_bis(2,:,1)=u(tri(deformed(i),:));
        uv_tricorner_bis(1,:,2)=v(tri(deformed(i),:));
        uv_tricorner_bis(2,:,2)=v(tri(deformed(i),:));
        
        xy_tricorner_bis(1,j_other,1)=x_test;
        xy_tricorner_bis(2,j_other_end,1)=x_test;
        xy_tricorner_bis(1,j_other,2)=y_test;
        xy_tricorner_bis(2,j_other_end,2)=y_test;
        
        uv_tricorner_bis(1,j_other,1)=u_test;
        uv_tricorner_bis(2,j_other_end,1)=u_test;
        uv_tricorner_bis(1,j_other,2)=v_test;
        uv_tricorner_bis(2,j_other_end,2)=v_test;
        
        [data_bis.dudx,data_bis.dudy,data_bis.dvdx,data_bis.dvdy,data_bis.scale,data_bis.area]=partial_deriv(xy_tricorner_bis,uv_tricorner_bis);
        [invar_bis]=invariants(data_bis.dudx,data_bis.dudy,data_bis.dvdx,data_bis.dvdy);
        
        if(max(abs(invar_bis.div))<=abs(div_min_max))
            div_min_max=max(abs(invar_bis.div));
            
            j_other_best=j_other;
            j_other_end_best=j_other_end;
            
            x_best=x_test;
            y_best=y_test;
            u_best=u_test;
            v_best=v_test;
        end
        
    end
    
    tri_new(deformed(i),j_other_best)=nb_init+i;
    tri_new(nb_init+i,:)=tri(deformed(i),:);
    tri_new(nb_init+i,j_other_end_best)=nb_init+i;
    
    x_new(nb_init+i)=x_best;
    y_new(nb_init+i)=y_best;
    u_new(nb_init+i)=u_best;
    v_new(nb_init+i)=v_best;
    
    %
    %
    %
    %     end
    %
    %     if(j_best>0)
    %         tri_new(deformed(i),j_other_best)=nb_init+i;
    %         tri_new(nb_init+i,:)=tri(deformed(i),:);
    %         tri_new(nb_init+i,j_best)=nb_init+i;
    %
    %         x_new(nb_init+i)=x_best;
    %         y_new(nb_init+i)=y_best;
    %         u_new(nb_init+i)=u_best;
    %         v_new(nb_init+i)=v_best;
    %     end
end
tri_new=tri_new(tri_new(:,1)>0,:);

% deformation with the triangulation
xy_new_tricorner(:,:,1)=x_new(tri_new);
xy_new_tricorner(:,:,2)=y_new(tri_new);
uv_new_tricorner(:,:,1)=u_new(tri_new);
uv_new_tricorner(:,:,2)=v_new(tri_new);

[data_new.dudx,data_new.dudy,data_new.dvdx,data_new.dvdy,data_new.scale,data_new.area]=partial_deriv(xy_new_tricorner,uv_new_tricorner);

[invar_new]=invariants(data_new.dudx,data_new.dudy,data_new.dvdx,data_new.dvdy);

[k]=plots(tri_new,invar_new,angle_i,xc,yc,x_new,u_new,y_new,v_new,dx,dy,k,[],[]);
%
% figure
% mean(invar_tris.div)
% mean(abs(invar_tris.div))
% hist(invar_tris.div,[-0.03:0.0001:0.03])
% figure
% mean(invar_tris.shear)
% hist(invar_tris.shear,[0:0.0001:0.08])
% mean(invar_tris.eps)


end

function [result]=recursive_detect_feature(visited_triangle,triangle,level,tr,is_deformed,max_level)
next_triangle=neighbors(tr,triangle);
if(sum(isnan(next_triangle))==3) % a triangle alone in the blue (not taken into account)
    disp(tr)
    disp(next_triangle)
    disp(triangle)
    error('still a triangle alone in the blue !!!');
else
    visited_triangle=[visited_triangle,triangle];
    result=triangle;
end
if(level>max_level )
    return;
else
    for k=1:3,
        next_ti=next_triangle(k);
        if(~isnan(next_ti) && sum(next_ti==visited_triangle)==0)
            if(is_deformed(next_ti))
                [tmp_result]=recursive_detect_feature(visited_triangle,next_ti,level+1,tr,is_deformed,max_level);
                result=[result,tmp_result];
            end
        end
    end
end
end

% function [result,visited]=recursive_detect_feature(visited_triangle,triangle,level,tr,visited,max_level)
% result=triangle;
% visited_triangle=[visited_triangle,triangle];
% %visited(triangle)=1;
% if(level>max_level )
%     return;
% else
%     next_triangle=neighbors(tr,triangle);
%     for k=1:3,
%         next_ti=next_triangle(k);
%         if(~isnan(next_ti) && sum(next_ti==visited_triangle)==0)
%             if(~visited(next_ti))
%                 [tmp_result,visited]=recursive_detect_feature(visited_triangle,next_ti,level+1,tr,visited,max_level);
%                 result=[result,tmp_result];
%             end
%         end
%     end
% end
% end

function [error_div,error_conv,k,ex_deformed]=smart_smoothing(showplot,vec_max_level,data,invar,tr,x_center,y_center,tri,angle_i,xc,yc,x,u,y,v,dx,dy,k,scaling_displ,length_fract,min_def,ex_deformed)

%------------ second test correction --------------
deformed_bis=find(invar.eps>1e-9);

for ind_max_level=1:length(vec_max_level)
    max_level=vec_max_level(ind_max_level);
    
    data_tris=data;
    data_tris.size_kernel=zeros(size(data_tris.dudx));
   
    def=invar.eps';
    deformed=find(def>=min_def);
    
    is_deformed=zeros(size(invar.eps));
    is_deformed(deformed)=1;
    
    [min_dist,ind_min]=min(hypot(x_center(deformed)-0.5,y_center(deformed)-0.5));
    if(ex_deformed==-1)
        ex_deformed=deformed(ind_min);
    end
    
    for i=1:length(deformed)
        
        other_deformed=recursive_detect_feature([],deformed(i),1,tr,is_deformed,max_level);
        other_deformed=unique(other_deformed);
        if(deformed(i)==ex_deformed)
            ex_other_deformed=other_deformed;
        end
        
        dudx=0;
        dudy=0;
        dvdx=0;
        dvdy=0;
        for j=1:length(other_deformed)
            dudx=dudx+data.dudx(other_deformed(j))*data.area(other_deformed(j));
            dudy=dudy+data.dudy(other_deformed(j))*data.area(other_deformed(j));
            dvdx=dvdx+data.dvdx(other_deformed(j))*data.area(other_deformed(j));
            dvdy=dvdy+data.dvdy(other_deformed(j))*data.area(other_deformed(j));
        end
        
        dudx=dudx/sum(data.area(other_deformed(:)));
        dudy=dudy/sum(data.area(other_deformed(:)));
        dvdx=dvdx/sum(data.area(other_deformed(:)));
        dvdy=dvdy/sum(data.area(other_deformed(:)));
        
        
        theta=[0:0.01:pi];
        l=x_center(other_deformed)*cos(theta)+y_center(other_deformed)*sin(theta);
        std_l=std(l,0,1);
        [std_min,ind_i]=min(std_l);
        slope(i)=theta(ind_i);
        l_origin(i)=mean(l(:,ind_i));
        
        data_tris.dudx(deformed(i))=dudx;
        data_tris.dudy(deformed(i))=dudy;
        data_tris.dvdx(deformed(i))=dvdx;
        data_tris.dvdy(deformed(i))=dvdy;
        
        data_tris.size_kernel(deformed(i))=length(other_deformed);
    end
    
%     max(data_tris.size_kernel)
%     min(data_tris.size_kernel)
    disp(max_level)
    max(data_tris.size_kernel)/max_level^2
    
    [invar_tris]=invariants(data_tris.dudx,data_tris.dudy,data_tris.dvdx,data_tris.dvdy);
    
    figure(111)
    hold on
    %plot(max_level,max(invar_tris.shear),'.')
    plot(max_level,1/(sum(invar_tris.shear(deformed_bis).*data.area(deformed_bis))/sum(data.area(deformed_bis))),'.')
    
    if(length(dx)==1)
        error_div(ind_max_level)=sum(invar_tris.div(invar_tris.div>0).*data.area(invar_tris.div>0)) -length_fract(1)*dy(1)*(dy(1)>0);
        error_conv(ind_max_level)=sum(invar_tris.div(invar_tris.div<0).*data.area(invar_tris.div<0))-length_fract(1)*dy(1)*(dy(1)<0);
    else
        error_div(ind_max_level)=sum(invar_tris.div(invar_tris.div>0).*data.area(invar_tris.div>0)) -(length_fract(1)*dy(1)*(dy(1)>0)+length_fract(2)/2*dx(2)*(dx(2)>0)*sin(angle_i(2)-angle_i(1)));
        error_conv(ind_max_level)=sum(invar_tris.div(invar_tris.div<0).*data.area(invar_tris.div<0))-(length_fract(1)*dy(1)*(dy(1)<0)+length_fract(2)/2*dx(2)*(dx(2)<0)*sin(angle_i(2)-angle_i(1)));
    end
    
end


if(showplot)
    [k]=plots(tri,invar_tris,angle_i,xc,yc,x,u,y,v,dx,dy,k,scaling_displ,ex_deformed,ex_other_deformed);
end
end

function [k]=plots(tri,invar,angle_i,xc,yc,x,u,y,v,dx,dy,k,scaling_displ,ex_deformed,other_deformed)

figure
%subaxis(2,2,k, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05);
f=patch(x(tri)'+scaling_displ*u(tri)',y(tri)'+scaling_displ*v(tri)',invar.div','FaceColor','flat');
hold on
f=patch(x(tri(other_deformed,:))'+scaling_displ*u(tri(other_deformed,:))',y(tri(other_deformed,:))'+scaling_displ*v(tri(other_deformed,:))',invar.div(other_deformed,:)','FaceColor','flat','EdgeColor','w','LineWidth',2);
f=patch(x(tri(ex_deformed,:))'+scaling_displ*u(tri(ex_deformed,:))',y(tri(ex_deformed,:))'+scaling_displ*v(tri(ex_deformed,:))',invar.div(ex_deformed,:)','FaceColor','flat','EdgeColor','g','LineWidth',2);

xtest=[0:0.001:1];
ytest=min(1,max(0,yc(1)+tan(angle_i(1))*([0:0.001:1]-xc(1))));
ind=find((ytest>0).*(ytest<1));
plot(xtest(ind),ytest(ind),'.-')
ltest=[0:0.001:1];
if(length(dx)>1)
    xtest=min(1,xc(2)+cos(angle_i(2))*[0:0.001:1]);
    ytest=min(1,yc(2)+sin(angle_i(2))*[0:0.001:1]);
    ind=find((xtest>0).*(xtest<1).*(ytest>0).*(ytest<1));
    plot(xtest(ind),ytest(ind),'.-')
end

colorbar
%set(f,'edgecolor','none','facecolor','flat')
caxis([-0.04 0.04])
%caxis([-10*max(abs([dx,dy])),10*max(abs([dx,dy]))])
axis([0-max(abs([dx,dy]))*scaling_displ*2,1+max(abs([dx,dy]))*scaling_displ*2,0-max(abs([dx,dy]))*scaling_displ*2,1+max(abs([dx,dy]))*scaling_displ*2])
colormap(blue2red)
%cbfreeze; axis equal;
%freezeColors

figure
%subaxis(2,2,k+2, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.05);
f=patch(x(tri)'+scaling_displ*u(tri)',y(tri)'+scaling_displ*v(tri)',invar.shear','FaceColor','flat');
hold on
f=patch(x(tri(other_deformed,:))'+scaling_displ*u(tri(other_deformed,:))',y(tri(other_deformed,:))'+scaling_displ*v(tri(other_deformed,:))',invar.shear(other_deformed,:)','FaceColor','flat','EdgeColor','w','LineWidth',2);
f=patch(x(tri(ex_deformed,:))'+scaling_displ*u(tri(ex_deformed,:))',y(tri(ex_deformed,:))'+scaling_displ*v(tri(ex_deformed,:))',invar.shear(ex_deformed,:)','FaceColor','flat','EdgeColor','g','LineWidth',2);


xtest=[0:0.001:1];
ytest=min(1,max(0,yc(1)+tan(angle_i(1))*([0:0.001:1]-xc(1))));
ind=find((ytest>0).*(ytest<1));
plot(xtest(ind),ytest(ind),'.-')
ltest=[0:0.001:1];
if(length(dx)>1)
    xtest=min(1,xc(2)+cos(angle_i(2))*[0:0.001:1]);
    ytest=min(1,yc(2)+sin(angle_i(2))*[0:0.001:1]);
    ind=find((xtest>0).*(xtest<1).*(ytest>0).*(ytest<1));
    plot(xtest(ind),ytest(ind),'.-')
end
    
colorbar
%set(f,'edgecolor','none','facecolor','flat')
caxis([0 0.08])
%caxis([0,20*max(abs([dx,dy]))])
axis([0-max(abs([dx,dy]))*scaling_displ*2,1+max(abs([dx,dy]))*scaling_displ*2,0-max(abs([dx,dy]))*scaling_displ*2,1+max(abs([dx,dy]))*scaling_displ*2]);
colormap(gray2red)
%cbfreeze; axis equal;
%freezeColors

end

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

end

function [dudx,dudy,dvdx,dvdy]=sub_diff_contour(x,y,u,v,aire)
%function [dudx,dudy,dvdx,dvdy]=sub_diff_contour_x,y,u,v,aire)
%compute partial derivates (x et y coordinates have to be ordered in the counter-clockwise sense)

temp = 1./(2.*aire);

dudx=  temp.*( (u(:,2)+u(:,1)).*(y(:,2)-y(:,1)) + (u(:,3)+u(:,2)).*(y(:,3)-y(:,2)) + (u(:,1)+u(:,3)).*(y(:,1)-y(:,3)) );
dudy= -temp.*( (u(:,2)+u(:,1)).*(x(:,2)-x(:,1)) + (u(:,3)+u(:,2)).*(x(:,3)-x(:,2)) + (u(:,1)+u(:,3)).*(x(:,1)-x(:,3)) );

dvdx=  temp.*( (v(:,2)+v(:,1)).*(y(:,2)-y(:,1)) + (v(:,3)+v(:,2)).*(y(:,3)-y(:,2)) + (v(:,1)+v(:,3)).*(y(:,1)-y(:,3)) );
dvdy= -temp.*( (v(:,2)+v(:,1)).*(x(:,2)-x(:,1)) + (v(:,3)+v(:,2)).*(x(:,3)-x(:,2)) + (v(:,1)+v(:,3)).*(x(:,1)-x(:,3)) );

return

end
