function [out]=EGPS_drift_deformation(startdate_str,enddate_str,compute_deformation,limit_domain)
% Uses the EGPS image pairs dataset downloaded from the Ron Kwok website and calculates sea ice
% deformations from the image pairs.
%
% example: EGPS_drift_deformation('2008_03_01','2008_03_05',0,[])

% flags
test=0;

if(compute_deformation)
    % parameters to reject small datastream
    min_size_stream=200;    % typical value = 200       for RGPS (to be tested for EGPS)
    
    % parameters to reject bad triangles after the triangulation
    min_area=5*1e6;         % typical value = 5*1e6     for RGPS (To be tested for EGPS)
    max_area=400*1e6;       % typical value = 400*1e6   for RGPS (To be tested for EGPS)
    min_minang=5;           % typical value = 5         for RGPS (To be tested for EGPS)
    max_long_side=25*1e3;   % typical value = 25*1e3    for RGPS (To be tested for EGPS)

    % % parameters for the smart smoother
    min_def=0.02; % std value = 0.02
    max_level=3; % std value =3
    
    filter_noise=1;
    if(test==1)
        filter_noise=1
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

szdat = 1000000;

% Calculate derivatives and strain rates from lagriangian cells
% (triangles).
% Process each image pair seperatly.
% Calculation of area A and derivatives after:
% RADARSAT GEOPHYSICAL PROCESSOR SYSTEM, DATA USER???S HANDBOOK (Version 1.0)
% or Kwok et al. (2008), Variability of sea ice simulations assessed with RGPS kinematics, JGR, doi:10.1029/2008JC004783

out.xy = zeros(2*szdat,2);
out.uv = zeros(2*szdat,2);
out.dnum = zeros(2*szdat,1);
out.deltat = zeros(2*szdat,1);
out.stream = char(ones(2*szdat,1)*'z');     % no stream in EGPS, set stream equals to 'z'

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
        outfile = ['EGPS_deformation_min_' num2str(min_size_stream) '_eps_00' num2str(min_def*100) '_level_' num2str(max_level) '.mat'];
    else
        outfile = ['EGPS_deformation_min_' num2str(min_size_stream) '_no_filter.mat'];
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

files=[];

if(startyear==endyear)
    tmp_year=num2str(startyear);
    for month=str2num(startmonth):str2num(endmonth)
        folder=dir([indir 'EGPS_ice_drift/EULIMGPAIR/cenarc_' , tmp_year(3:4) , sprintf('%02d',month), ' *']);
        files=[dir([indir 'EGPS_ice_drift/EULIMGPAIR/', folder.name, '/*.hdf']); files];
    end
else
    tmp_year=num2str(startyear);
    for month=str2num(startmonth):12
        folder=dir([indir 'EGPS_ice_drift/EULIMGPAIR/cenarc_' , tmp_year(3:4) , sprintf('%02d',month), ' *']);
        files=[dir([indir 'EGPS_ice_drift/EULIMGPAIR/', folder.name, '/*.hdf']); files];
    end
    tmp_year=num2str(endyear);
    for month=1:str2num(endmonth)
        folder=dir([indir 'EGPS_ice_drift/EULIMGPAIR/cenarc_' , tmp_year(3:4) , sprintf('%02d',month), ' *']);
        files=[dir([indir 'EGPS_ice_drift/EULIMGPAIR/', folder.name, '/*.hdf']); files];
    end
end

first_image_pair=1;
last_image_pair=length(files);

if(last_image_pair == 0)
    disp('No input files found!')
    return
end

jend=0;
for i=first_image_pair:last_image_pair
    
    file=files(i);
    filename=file.name;
    
    date_before=datenum(filename(1:12),'yyyymmddHHMMSS');
    date=       datenum(filename(14:25),'yyyymmddHHMMSS');
    
    deltat = date - date_before; % day
    dnum = date;
    
    if(date_before>dnum_start && date<dnum_end && date_before<date)
        
        FINFO=hdfread(filename,'NGRID With Obs:');
        n2=str2num(FINFO{1}');
       
        lat=hdfread(filename,'src_grid_lat(deg)')';
        lon=hdfread(filename,'src_grid_lon(deg)')';
        
        m_proj('Stereographic','lon',-45,'lat',90,'radius',60);
        [x,y]=m_ll2xy(lon,lat);
        x=(x*6378.273)*1000; % m
        y=(y*6378.273)*1000; % m
        
        trux2=hdfread(filename,'grid_dx(km)')*1000/deltat; % m/day
        trvy2=hdfread(filename,'grid_dy(km)')*1000/deltat; % m/day
        
        % Limit the domain
        if(length(limit_domain)==4)
            ind_selected=find((x>=limit_domain(1)).*(x<=limit_domain(2)).*(y>=limit_domain(3)).*(y<=limit_domain(4)));
            if(length(ind_selected)==0)
                continue;
            end
            x=x(ind_selected);
            y=y(ind_selected);
            trux2=trux2(ind_selected);
            trvy2=trvy2(ind_selected);
        end
                
        % treatement of the data
        txy2=[x,y]; % m
        n_selected=length(trux2);
        
        tdeltat2=deltat;
        tdnum2=dnum;
        
        if(compute_deformation)
        if n2>min_size_stream
                 
            [ntri,ttdnum,ttdeltat,area,xy_tricorner,uv,uv_tricorner,dudx,dudy,dvdx,dvdy,id_stream,quality_index]=def_from_drift(txy2,tdeltat2,tdnum2,trux2,trvy2,'z',min_area,max_area,min_minang,max_long_side,showplot,filter_noise,min_def,max_level);
            
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
            % add to output structure
            jend = jstart+n_selected-1;
            
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
