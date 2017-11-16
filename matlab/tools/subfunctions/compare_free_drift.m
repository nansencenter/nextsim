function [u1_select,...
    v1_select,...
    speed1_select,...
    x1_select,...
    y1_select,...
    u2_select,...
    v2_select,...
    speed2_select,...
    x2_select,...
    y2_select,...
    u3_select,...
    v3_select,...
    speed3_select,...
    x3_select,...
    y3_select,...
    u4_select,...
    v4_select,...
    speed4_select,...
    x4_select,...
    y4_select,...
    i_box_select,...
    j_box_select,...
    xmin,...
    ymin,...
    min_box_size,...
    ratio_dom_resol]=compare_free_drift(defo_vec,title_defo)
% example:
% compare_free_drift({'defo_mobs_test12_rgps.mat','defo_mobs_test12_freedrift_rgps.mat','defo_rgps.mat'},{'Model','obs'});

u1_select=[];
v1_select=[];
speed1_select=[];
x1_select=[];
y1_select=[];
u2_select=[];
v2_select=[];
speed2_select=[];
x2_select=[];
y2_select=[];
u3_select=[];
v3_select=[];
speed3_select=[];
x3_select=[];
y3_select=[];
u4_select=[];
v4_select=[];
speed4_select=[];
x4_select=[];
y4_select=[];

if(length(defo_vec)~=3 && length(defo_vec)~=4)
    error('The length of the input vector should be 3 or 4')
end

% Uses the same range as in MacPhee (1979) 8cm/s to 22 cm/s
min_speed= 8/100/1000*3600*24; % km/day
max_speed=22/100/1000*3600*24; % km/day

% parameters that could be changed
plot_the_consistency_plots=0;   % set to 1 if you want to check that the temporal and spatial mapping are right

% load of the data
[defo_1,ebox_1]=load_defo(defo_vec{1});
[defo_2,ebox_2]=load_defo(defo_vec{2});   
[defo_3,ebox_3]=load_defo(defo_vec{3}); 
if(length(defo_vec)==4)
    [defo_4,ebox_4]=load_defo(defo_vec{4});    
end

if(isempty(ebox_1) || isempty(ebox_2)||isempty(ebox_3))
    return;
end
if(length(defo_vec)==4 && isempty(ebox_4))
    return;
end

data_1=defo_1.data;
data_2=defo_2.data;
data_3=defo_3.data;
if(length(defo_vec)==4)
    data_4=defo_4.data;
end

xmin=ebox_1.xmin;
ymin=ebox_1.ymin;
min_box_size=ebox_1.min_box_size;
ratio_dom_resol=ebox_1.ratio_dom_resol;

% selection of the boxes (here we take all the domain)
tmp_ebox_1=ebox_1.mask;
tmp_ebox_2=ebox_2.mask;
tmp_ebox_3=ebox_3.mask;
if(length(defo_vec)==4)
    tmp_ebox_4=ebox_4.mask;
end

% search the corresponding indices
if(length(defo_vec)==4)
    [i_tmp,j_tmp]=find(tmp_ebox_1.*tmp_ebox_2.*tmp_ebox_3.*tmp_ebox_4);
else
    [i_tmp,j_tmp]=find(tmp_ebox_1.*tmp_ebox_2.*tmp_ebox_3);
end

% Loop over the data
indices_1=[];
indices_2=[];
indices_3=[];
if(length(defo_vec)==4)
    indices_4=[];
end
i_box=[];
j_box=[];

for i=1:length(i_tmp)

    tmp_indices_1=ebox_1.full{i_tmp(i),j_tmp(i)}';
    out_1_stream=char(data_1.stream(tmp_indices_1));
    out_1_dnum=data_1.dnum(tmp_indices_1);
    out_1_deltat=data_1.deltat(tmp_indices_1);
    
    tmp_indices_2=ebox_2.full{i_tmp(i),j_tmp(i)}';
    out_2_stream=char(data_2.stream(tmp_indices_2));
    out_2_dnum=data_2.dnum(tmp_indices_2);
    out_2_deltat=data_2.deltat(tmp_indices_2);
    
    tmp_indices_3=ebox_3.full{i_tmp(i),j_tmp(i)}';
    out_3_stream=char(data_3.stream(tmp_indices_3));
    out_3_dnum=data_3.dnum(tmp_indices_3);
    out_3_deltat=data_3.deltat(tmp_indices_3);
    
    if(length(defo_vec)==4)
        tmp_indices_4=ebox_4.full{i_tmp(i),j_tmp(i)}';
        out_4_stream=char(data_4.stream(tmp_indices_4));
        out_4_dnum=data_4.dnum(tmp_indices_4);
        out_4_deltat=data_4.deltat(tmp_indices_4);
    end
    
    % for each stream
    streams=unique(out_1_stream);
    %for s=1:length(streams)
%         ind_1_k=find(out_1_stream==streams(s))
%         ind_2_k=find(out_2_stream==streams(s))
%         ind_3_k=find(out_3_stream==streams(s))
        
        ind_1_k=find(out_1_stream);
        ind_2_k=find(out_2_stream);
        ind_3_k=find(out_3_stream);
        if(length(defo_vec)==4)
            ind_4_k=find(out_4_stream);
        end

        if(length(defo_vec)==4)
            if(isempty(ind_2_k)||isempty(ind_3_k)||isempty(ind_4_k)) break; end
        else
            if(isempty(ind_2_k)||isempty(ind_3_k)) break; end
        end
        
        % for each date
        dnums=unique(out_1_dnum(ind_1_k));
        for l=1:length(dnums)
%             datestr(out_1_dnum(ind_1_k))
%             datestr(out_2_dnum(ind_2_k))
%             datestr(out_3_dnum(ind_3_k))
%             datestr(dnums(l))
            ind_1_l=find(out_1_dnum(ind_1_k)==dnums(l));
            ind_2_l=find(out_2_dnum(ind_2_k)==dnums(l));
            ind_3_l=find(out_3_dnum(ind_3_k)==dnums(l));
            if(length(defo_vec)==4)
                ind_4_l=find(out_4_dnum(ind_4_k)==dnums(l));
            end
            
            if(length(defo_vec)==4)
                if(isempty(ind_2_l)||isempty(ind_3_l)||isempty(ind_4_l)) break; end
            else
                if(isempty(ind_2_l)||isempty(ind_3_l)) break; end
            end
            
            % for each deltat
            deltats=unique(out_1_deltat(ind_1_k(ind_1_l)));
            for m=1:length(deltats)
                ind_1_m=find(out_1_deltat(ind_1_k(ind_1_l))==deltats(m));
                ind_2_m=find(out_2_deltat(ind_2_k(ind_2_l))==deltats(m));
                ind_3_m=find(out_3_deltat(ind_3_k(ind_3_l))==deltats(m));
                if(length(defo_vec)==4)
                    ind_4_m=find(out_4_deltat(ind_4_k(ind_4_l))==deltats(m));
                end
                
                if(length(defo_vec)==4)
                    if(isempty(ind_2_m)||isempty(ind_3_m)||isempty(ind_4_m)) break; end
                else
                    if(isempty(ind_2_m)||isempty(ind_3_m)) break; end
                end
                
                if(length(defo_vec)==4)
                    nb_match=min(min(min(length(ind_1_m),length(ind_2_m)),length(ind_3_m)),length(ind_4_m));
                else
                    nb_match=min(min(length(ind_1_m),length(ind_2_m)),length(ind_3_m));
                end
                tmp_indices_1(ind_1_k(ind_1_l(ind_1_m(1:nb_match))));
                indices_1=[indices_1;tmp_indices_1(ind_1_k(ind_1_l(ind_1_m(1:nb_match))))];
                indices_2=[indices_2;tmp_indices_2(ind_2_k(ind_2_l(ind_2_m(1:nb_match))))];
                indices_3=[indices_3;tmp_indices_3(ind_3_k(ind_3_l(ind_3_m(1:nb_match))))];
                if(length(defo_vec)==4)
                    indices_4=[indices_4;tmp_indices_4(ind_4_k(ind_4_l(ind_4_m(1:nb_match))))];  
                end
                i_box=[i_box,i_tmp(i)*ones(1,nb_match)];
                j_box=[j_box,j_tmp(i)*ones(1,nb_match)];
            end
        end
%    end
end

if(isempty(indices_1))
    disp('empty')
    return
end
    
% get the data
[u1,v1,x1,y1,dnum1,speed1]=get_data(data_1,indices_1);
[u2,v2,x2,y2,dnum2,speed2]=get_data(data_2,indices_2);
[u3,v3,x3,y3,dnum3,speed3]=get_data(data_3,indices_3);
if(length(defo_vec)==4)
    [u4,v4,x4,y4,dnum4,speed4]=get_data(data_4,indices_4);
end

if(exist('speed4','var'))
    ind=find((speed3>=min_speed).*(speed3<=max_speed).*(speed4>=min_speed).*(speed4<=max_speed).*(speed1>=min_speed).*(speed1<=max_speed));
else
    ind=find((speed3>=min_speed).*(speed3<=max_speed).*(speed1>=min_speed).*(speed1<=max_speed));
end

u1=u1(ind);
v1=v1(ind);
x1=x1(ind);
y1=y1(ind);
dnum1=dnum1(ind);
speed1=speed1(ind);
u2=u2(ind);
v2=v2(ind);
x2=x2(ind);
y2=y2(ind);
dnum2=dnum2(ind);
speed2=speed2(ind);
u3=u3(ind);
v3=v3(ind);
x3=x3(ind);
y3=y3(ind);
dnum3=dnum3(ind);
speed3=speed3(ind);
if(length(defo_vec)==4)
    u4=u4(ind);
    v4=v4(ind);
    x4=x4(ind);
    y4=y4(ind);
    dnum4=dnum4(ind);
    speed4=speed4(ind);
end

i_box=i_box(ind);
j_box=j_box(ind);

% copy dnum from the observations
dnum1=dnum3;
dnum2=dnum3;
dnum4=dnum3;

% plotting consistency plots
if(plot_the_consistency_plots)
    figure;
    plot(y1,y2,'.')
    
    figure;
    plot(x1,x2,'.')
    
    figure;
    plot(x1,y1,'.')
    
    figure;
    plot(dnum1,dnum2,'.')
    
    figure;
    plot(double(data_1.stream(indices_1,1)),double(data_2.stream(indices_2,1)),'.')
end

dif_ref_freedrift=sqrt((u2-u1).^2+(v2-v1).^2);
%rel_dif_ref_freedrift=dif_ref_freedrift./speed2;
%rel_dif_ref_freedrift=dif_ref_freedrift./speed1;
rel_dif_ref_freedrift=dif_ref_freedrift./speed3;
max_rel_dif=0.10;
ind_freedrift=find(rel_dif_ref_freedrift<max_rel_dif);

% % plot of the error in terms of velocity
% figure
% scatter(x1,y1,4,dif_ref_freedrift)
% colorbar
% caxis([0 max_speed_to_plot])
% 
% % plot of the error in terms of velocity
% figure
% scatter(x1,y1,4,rel_dif_ref_freedrift)
% colorbar
% caxis([0 max_rel_dif])

u1_select=u1(ind_freedrift);
v1_select=v1(ind_freedrift);
speed1_select=speed1(ind_freedrift);
x1_select=x1(ind_freedrift);
y1_select=y1(ind_freedrift);

u2_select=u2(ind_freedrift);
v2_select=v2(ind_freedrift);
speed2_select=speed2(ind_freedrift);
x2_select=x2(ind_freedrift);
y2_select=y2(ind_freedrift);

u3_select=u3(ind_freedrift);
v3_select=v3(ind_freedrift);
speed3_select=speed3(ind_freedrift);
x3_select=x3(ind_freedrift);
y3_select=y3(ind_freedrift);

if(length(defo_vec)==4)
    u4_select=u4(ind_freedrift);
    v4_select=v4(ind_freedrift);
    speed4_select=speed4(ind_freedrift);
    x4_select=x4(ind_freedrift);
    y4_select=y4(ind_freedrift);
end

i_box_select=i_box(ind_freedrift);
j_box_select=j_box(ind_freedrift);

end

function [u,v,x,y,dnum,speed]=get_data(data,indices);
u=data.uv(indices,1)/1000; % in km/day
v=data.uv(indices,2)/1000; % in km/day
x=data.x(indices)/1000;
y=data.y(indices)/1000;
dnum=data.dnum(indices,1);
speed=sqrt(u.^2+v.^2);

end