function [out]=DTU_drift_deformation_txt(startdate_str,enddate_str,compute_deformation,limit_domain)
% Uses the DTU sarbased dataset downloaded from the MyOcean website and calculates sea ice
% deformations from the drift vectors.
%
% example: DTU_drift_deformation('2015_10_01','2015_10_05',1,[])

% flags
test=1;

%folder_name='DTU_ice_drift';
folder_name='DTU_ice_drift_100m';

if(compute_deformation)
    % parameters to reject small datastream
    min_size_stream=200;
    
    % parameters to reject bad triangles after the triangulation
    min_area=5*1e6;         % typical value = 5*1e6     for RGPS (also ok for Globice)
    max_area=400*1e6;       % typical value = 400*1e6   for RGPS (also ok for Globice)
    min_minang=5;           % typical value = 5         for RGPS (also ok for Globice)
    max_long_side=25*1e3;   % typical value = 25*1e3    for RGPS (also ok for Globice)

    % % parameters for the smart smoother
    min_def=0.02; % std value = 0.02
    max_level=3; % std value =3
    
    filter_noise=1;
    if(test==1)
        filter_noise=0
    end
end

define_default_data_path;

startdate=datestr(datenum(startdate_str),'yyyy_mm_dd')
enddate  =datestr(datenum(  enddate_str),'yyyy_mm_dd')

startyear = startdate(1:4);
startmonth = startdate(6:7);
startday = startdate(9:10);
endyear = enddate(1:4);
endmonth = enddate(6:7);
endday = enddate(9:10);

dnum_start=datenum(startdate);
dnum_end=datenum(enddate);

outdir = indir;

showplot = test;    % show some extra triangulation plots (slows down program)

% ---- End user variables --------------

% Initialization
szdat = 1000000;

out.xy = zeros(2*szdat,2);
out.uv = zeros(2*szdat,2);
out.dnum = zeros(2*szdat,1);
out.deltat = zeros(2*szdat,1);
out.stream = char(ones(2*szdat,1)*'z');     % no stream in Globice, set stream equals to 'z'

if(compute_deformation)
    out.dudx = zeros(2*szdat,1);
    out.dudy = zeros(2*szdat,1);
    out.dvdx = zeros(2*szdat,1);
    out.dvdy = zeros(2*szdat,1);
    out.area = zeros(2*szdat,1);
    out.xy_tricorner = zeros(2*szdat,3,2);
    out.uv_tricorner = zeros(2*szdat,3,2);
    out.size_datatake = zeros(2*szdat,1);
    out.quality_index_tot=zeros(2*szdat,1);
    
    if(filter_noise)
        outfile = ['DTU_deformation_min_' num2str(min_size_stream) '_eps_00' num2str(min_def*100) '_level_' num2str(max_level) '.mat'];
    else
        outfile = ['DTU_deformation_min_' num2str(min_size_stream) '_no_filter.mat'];
    end
end

% start loop
% Calculate derivatives and strain rates from the pairs of images.
% Process each image pair seperatly.
% Calculation of area A and derivatives after:
% RADARSAT GEOPHYSICAL PROCESSOR SYSTEM, DATA USER???S HANDBOOK (Version 1.0)
% or Kwok et al. (2008), Variability of sea ice simulations assessed with RGPS kinematics, JGR, doi:10.1029/2008JC004783
jstart = 1;
jjstart =1;

% Change of satelite at this pivot date
pivot_date=datenum('2012_05_01');
if((datenum(startdate_str)<=pivot_date).*(datenum(  enddate_str)<=pivot_date))
    satelite='wsm_envisat';
    length_name=39;
elseif((datenum(startdate_str)>=pivot_date).*(datenum(  enddate_str)>=pivot_date))
    satelite='sarbased';
    length_name=36;
else
    error(['There is a change of satellite on the ', datestr(pivot_date),'. Please chose a time interval taht do not contain this pivot date.']);
end

if(strcmp(folder_name,'DTU_ice_drift'))
    prefix_file=['ice_drift_polstereo_', satelite, '_north_']
else
    prefix_file='';
    length_name=1;
end
    
files=[];
if(startyear==endyear)
    for month=str2num(startmonth):str2num(endmonth)
        str_month=sprintf('%02d',month);
        files=[dir([indir, folder_name, '/', startyear, '/' prefix_file, startyear , str_month, '*']); files];
    end
else
    for month=str2num(startmonth):12
        str_month=sprintf('%02d',month);
        files=[dir([indir, folder_name '/', startyear, '/' prefix_file , startyear , str_month, '*']); files];
    end
    for month=1:str2num(endmonth)
        str_month=sprintf('%02d',month);
        files=[dir([indir folder_name '/', endyear,'/' prefix_file ,   endyear , str_month, '*']); files];
    end
end
    
first_image_pair=1;
last_image_pair=length(files);

if(last_image_pair == 0)
    [indir folder_name '/', endyear,'/*']
    disp('No input files found!')
    return
end

jend=0;
for i=first_image_pair:last_image_pair
    
    file=files(i);
    filename=file.name;

    if(strcmp(folder_name,'DTU_ice_drift'))
        date_before=datenum(filename(length_name   :length_name+13),'yyyymmddHHMMSS');
        date=       datenum(filename(length_name+15:length_name+28),'yyyymmddHHMMSS');
    else
        date_before=datenum(filename(1   :1+14),'yyyymmdd_HHMMSS');
        date=       datenum(filename(1+16:1+30),'yyyymmdd_HHMMSS');
    end    
        
    deltat = date - date_before; % day
    dnum = date;
    
    if(date_before>dnum_start && date<dnum_end && date_before<date)
        
        [dat] = read_DTUdrift(filename);
        
        % Limit the domain
        if(length(limit_domain)==4)
            ind_selected=find((dat.x*1000>=limit_domain(1)).*(dat.x*1000<=limit_domain(2)).*(dat.y*1000>=limit_domain(3)).*(dat.y*1000<=limit_domain(4)));
            if(length(ind_selected)==0)
                continue;
            end
            dat.dnum=dat.dnum(ind_selected);
            dat.x=dat.x(ind_selected);
            dat.y=dat.y(ind_selected);
            dat.rux=dat.rux(ind_selected);
            dat.rvy=dat.rvy(ind_selected);
            dat.deltat=dat.deltat(ind_selected);
        end
        
        tdnum2 = dat.dnum(1);
        txy2 = [dat.x, dat.y]*1000; % m
        trux2 = dat.rux*86400;    %m/day
        trvy2 = dat.rvy*86400;    %m/day
        tdeltat2 = dat.deltat(1);
         
        % Treatment of the data
        n2 = length(dat.dnum);
        if(compute_deformation)
            if n2>min_size_stream
  
                [ntri,ttdnum,ttdeltat,area,xy_tricorner,uv,uv_tricorner,dudx,dudy,dvdx,dvdy,id_stream,quality_index]=def_from_drift(txy2,tdeltat2,tdnum2,trux2,trvy2,'z',min_area,max_area,min_minang,max_long_side,showplot,filter_noise,min_def,max_level);
                
%                 if(showplot)
%                     figure(10000)
%                     scatter(x,y,36,qflag_1,'fill')
%                     figure(10001)
%                     scatter(x,y,36,qflag_2,'fill')
%                 end
                
                if(ntri>0)
                    % add to output structure
                    jend = jstart+ntri-1;
                    
                    tmp_length=length(out.dnum);
                    if(jend>tmp_length)
                        warning(['need to increase the size of out strcuture. You should increase szdat to at least ', num2str(jend)])
                        out.dnum=[out.dnum;zeros(tmp_length,1)];
                        out.deltat=[out.deltat;zeros(tmp_length,1)];
                        out.area=[out.area;zeros(tmp_length,1)];
                        out.xy_tricorner=[out.xy_tricorner;zeros(tmp_length,3,2)];
                        out.uv_tricorner=[out.uv_tricorner;zeros(tmp_length,3,2)];
                        out.uv=[out.uv;zeros(tmp_length,2)];
                        out.dudx=[out.dudx;zeros(tmp_length,1)];
                        out.dudy=[out.dudy;zeros(tmp_length,1)];
                        out.dvdx=[out.dvdx;zeros(tmp_length,1)];
                        out.dvdy=[out.dvdy;zeros(tmp_length,1)];
                        out.stream = [out.stream;char(ones(tmp_length,1)*'z')];
                        out.size_datatake=[out.size_datatake;zeros(tmp_length,1)];
                    end
                    
                    out.dnum(jstart:jend) = ttdnum;
                    out.deltat(jstart:jend) = ttdeltat;
                    out.area(jstart:jend) = area/1e6; % km^2
                    out.xy_tricorner(jstart:jend,:,:) = xy_tricorner/1000; % in km
                    out.uv(jstart:jend,:) = uv./1000; % in km/day
                    out.uv_tricorner(jstart:jend,:,:) = uv_tricorner/1000;
                    out.dudx(jstart:jend) = dudx; out.dudy(jstart:jend) = dudy; % in /day
                    out.dvdx(jstart:jend) = dvdx; out.dvdy(jstart:jend) = dvdy; % in /day
                    out.stream(jstart:jend) = id_stream;
                    out.size_datatake(jstart:jend) = ntri;
                    
                    if(jjstart>length(out.quality_index_tot))
                        out.quality_index_tot=[out.quality_index_tot;zeros(length(out.quality_index_tot),1)];
                    end
                    out.quality_index_tot(jjstart)=quality_index;
                    
                    jstart = jend+1;
                    jjstart=jjstart+1;
                end
            end
        else
            % velocities only
            
            % add to output structure
            jend = jstart+n2-1;
                
            tmp_length=length(out.dnum);
            if(jend>tmp_length)
                warning(['need to increase the size of out strcuture. You should increase szdat to at least ', num2str(jend)])
                out.dnum=[out.dnum;zeros(tmp_length,1)];
                out.deltat=[out.deltat;zeros(tmp_length,1)];
                out.xy=[out.xy;zeros(tmp_length,2)];
                out.uv=[out.uv;zeros(tmp_length,2)];
                out.stream = [out.stream;char(ones(tmp_length,1)*'z')];
            end
            
            out.dnum(jstart:jend) = tdnum2;
            out.deltat(jstart:jend) = tdeltat2;
            out.xy(jstart:jend,:) = txy2./1000; % in km
            out.uv(jstart:jend,1) = trux2./1000; % in km/day
            out.uv(jstart:jend,2) = trvy2./1000; % in km/day
            out.stream(jstart:jend) = 'z';

            jstart = jend+1;
            jjstart=jjstart+1;            
        end
    end
end % image pairs

% shorten dataset to actual length
out.dnum = out.dnum(1:jend);
out.deltat = out.deltat(1:jend);
out.xy = out.xy(1:jend,:);
out.uv = out.uv(1:jend,:);
out.stream = out.stream(1:jend);

if(compute_deformation)
    out.area = out.area(1:jend);
    out.dudx = out.dudx(1:jend);
    out.dudy = out.dudy(1:jend);
    out.dvdx = out.dvdx(1:jend);
    out.dvdy = out.dvdy(1:jend);
    out.xy_tricorner = out.xy_tricorner(1:jend,:,:);
    out.uv_tricorner = out.uv_tricorner(1:jend,:,:);
    out.quality_index_tot=out.quality_index_tot(1:jjstart-1,:);
    out.size_datatake = out.size_datatake(1:jend);

    % divergence
    out.rdiv = out.dudx + out.dvdy;
    % vorticity
    out.rvor = out.dvdx - out.dudy;
    % shear
    out.rshr = sqrt((out.dudx-out.dvdy).^2+(out.dudy+out.dvdx).^2);
    
    disp(['Save data to ' outfile])
    save(outfile,'out');
end

end
