function [data_bin,data_bin_mean]=coarsegrainer(ebox,data,nb_bins,min_box_size,min_area_factor)

if(isempty(ebox)||isempty(data)||isempty(data.indices))
    data_bin=[];
    data_bin_mean=[];
    return;
end
    
k=nb_bins+1;
% Loop over the bins
indices = data.indices;
area = data.area(indices);
data_bin(k).dudx=data.dudx(indices);
data_bin(k).dudy=data.dudy(indices);
data_bin(k).dvdx=data.dvdx(indices);
data_bin(k).dvdy=data.dvdy(indices);
data_bin(k).x=data.x(indices);
data_bin(k).y=data.y(indices);
data_bin(k).scale=sqrt(area);
data_bin(k).area=area;
data_bin(k).u=data.uv(indices,1)/1000;
data_bin(k).v=data.uv(indices,2)/1000;
data_bin(k).speed=sqrt(data_bin(k).u.^2+data_bin(k).v.^2);
data_bin(k).deltat=data.deltat(indices);

ebox_not_empty=ebox.mask;

% conc, damage, thickness only for the lowest scale
if(isfield(data,'damage'))
    if(~isempty(data.damage))
    data_bin(k).damage=data.damage(indices);
    data_bin(k).conc=data.conc(indices);
    data_bin(k).thick=data.thick(indices);
    end
end

data_u=data.uv(:,1)';
data_v=data.uv(:,2)';

tot_size=2^(nb_bins-1);

marsan_categories=[13,20,40,80,160,320,640,1280,2560,5120,10240];

for k=1:1:nb_bins
    % Loop over the boxes of each box size except the smallest one 10km for
    % the RGPS
    box_width = 2^(nb_bins-k);
    dist_boxcenter= 2^(nb_bins-1-k);
    
    min_area=(min_box_size*box_width)^2*min_area_factor*1e6;
    
    if(k==nb_bins)
        nb_box=tot_size;
    else
        nb_box=tot_size/dist_boxcenter;
    end
    
    l_max=nb_box*nb_box;
        
    deltat_tmp=zeros(l_max,1);
    dudx_tmp=zeros(l_max,1);
    dudy_tmp=zeros(l_max,1);
    dvdx_tmp=zeros(l_max,1);
    dvdy_tmp=zeros(l_max,1);
    x_tmp=zeros(l_max,1);
    y_tmp=zeros(l_max,1);
    area_tmp=zeros(l_max,1);
    u_tmp=zeros(l_max,1);
    v_tmp=zeros(l_max,1);
        
    % L counts the number of non-empty boxes
    L=0;
    
    for i=1:nb_box,
        for j=1:nb_box,

            if(k==nb_bins)
                look_i=i;
                look_j=j;
                
                not_empty=ebox_not_empty(look_i,look_j);
            else
                max_i=min((i-1)*dist_boxcenter+box_width,tot_size);
                max_j=min((j-1)*dist_boxcenter+box_width,tot_size);
            
                look_i=(i-1)*dist_boxcenter+1:max_i;
                look_j=(j-1)*dist_boxcenter+1:max_j;
                
                not_empty=sum(sum(ebox_not_empty(look_i,look_j)));
            end
                   
            if(not_empty>=box_width)            
                
                
                indices=[];
                for tmp_i=1:length(look_i),
                    for tmp_j=1:length(look_j)
                        indices=[indices,ebox.no_overlap{look_i(tmp_i),look_j(tmp_j)}];
                    end
                end
            
            % area, scale and derivative over the box
            % Empty boxes are not taken into account
            % keep boxes if the surface covered by RGPS cells or model elements is larger than
            % half the box area (min_area)
                area = sum(data.area(indices));
                if(area>min_area)
                    L=L+1;
                    n_element=length(indices);
                    weigth=data.area(indices)/area;
                    deltat_tmp(L,1) = sum(data.deltat(indices))/n_element;
                    dudx_tmp(L,1)=data.dudx(indices)'*weigth;
                    dudy_tmp(L,1)=data.dudy(indices)'*weigth;
                    dvdx_tmp(L,1)=data.dvdx(indices)'*weigth;
                    dvdy_tmp(L,1)=data.dvdy(indices)'*weigth;
                    x_tmp(L,1)=sum(data.x(indices))/n_element;
                    y_tmp(L,1)=sum(data.y(indices))/n_element;
                    area_tmp(L,1) = area;
                    u_tmp(L,1)=data_u(indices)*weigth/1000;
                    v_tmp(L,1)=data_v(indices)*weigth/1000;
                end
            end
        end
    end   
    data_bin(k).deltat=deltat_tmp(1:L,1); 
    data_bin(k).dudx=  dudx_tmp(1:L,1);
    data_bin(k).dudy=  dudy_tmp(1:L,1);
    data_bin(k).dvdx=  dvdx_tmp(1:L,1);
    data_bin(k).dvdy=  dvdy_tmp(1:L,1);
    data_bin(k).x=     x_tmp(1:L,1);
    data_bin(k).y=     y_tmp(1:L,1);
    data_bin(k).scale= sqrt(area_tmp(1:L,1));
    data_bin(k).area=  area_tmp(1:L,1); 
    data_bin(k).u=     u_tmp(1:L,1);
    data_bin(k).v=     v_tmp(1:L,1);
    data_bin(k).speed= sqrt(u_tmp(1:L,1).^2+v_tmp(1:L,1).^2);
end
    
q_vec=0.5*[0:6];
for k=1:(nb_bins+1)    
    % invariants of all the boxes
    data_bin(k).invar=invariants(data_bin(k).dudx,data_bin(k).dudy,data_bin(k).dvdx,data_bin(k).dvdy);     
    data_bin_mean.area(k)=mean(data_bin(k).area);
    data_bin_mean.length(k)=sqrt(data_bin_mean.area(k));
    
    
    for i=1:length(q_vec)
        q = q_vec(i);
        data_bin_mean.shear(k,i)=mean(data_bin(k).invar.shear.^q);
        data_bin_mean.div(k,i)=mean(abs(data_bin(k).invar.div).^q);
        data_bin_mean.eps(k,i)=mean(data_bin(k).invar.eps.^q);
        data_bin_mean.vor(k,i)=mean(abs(data_bin(k).invar.vor.^q));
    end
end

return
           