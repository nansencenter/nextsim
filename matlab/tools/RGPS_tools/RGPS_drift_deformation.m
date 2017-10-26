function [out]=RGPS_drift_deformation(startdate,enddate,compute_deformation,limit_domain)
% Uses the lagrangian RGPS sea ice motion
% dataset produced by RGPS_drift_vel_def.m and calculates sea ice
% deformations from the lagrangian cells.
%
% v1.0, Feb 2012, Gunnar Spreen
% v1.1, Mar 2012, added distance filter
% example: RGPS_drift_deformation('2006_12_03','2007_06_01',1,[])

% flags
test=0;

if(compute_deformation)
    % parameters to reject small datastream
    min_size_stream=200;

    % parameters to reject bad triangles after the triangulation
    min_area=5*1e6;
    max_area=400*1e6;
    min_minang=10;
    max_long_side=25*1e3;

    % % parameters for the smart smoother
    min_def=0.02; % std value = 0.02
    max_level=3; % std value =3
    
    filter_noise=1;
    if(test==1)
        filter_noise=1
    end
    
    if(~filter_noise)
        min_size_stream=3; % one triangle
    end
end
    
define_default_data_path;

files=dir([indir, 'RGPS_ice_drift/vel/RGPS_*']);
infile='';
for i=1:length(files)
    file_date_start=datenum(files(i).name( 6:15));
    file_date_end =datenum(files(i).name(17:26));
    if((file_date_start<=datenum(startdate)) && (file_date_end>=datenum(enddate)))
        infile = [indir 'RGPS_ice_drift/vel/' files(i).name];
    end
end

outdir = '';

showplot = test; % show some extra triangulation plots (slows down program)
quality = 0; % only use data with quality flag <=3 (not sure what that actually means)

% ---- End user variables --------------

if(exist(infile) ~= 2)
    disp(infile)
    disp('No input files found!')
    out.xy = zeros(1,2);
    out.uv = zeros(1,2);
    out.dnum = zeros(1,1);
    out.deltat = zeros(1,1);
    out.stream = char(ones(1,1)*'z');
    
    if(compute_deformation)
        out.dudx = zeros(1,1);
        out.dudy = zeros(1,1);
        out.dvdx = zeros(1,1);
        out.dvdy = zeros(1,1);
        out.area = zeros(1,1);
        out.xy_tricorner = zeros(1,3,2);
        out.uv_tricorner = zeros(1,3,2);
        out.size_datatake = zeros(1,1);
        out.quality_index_tot=zeros(1,1);
    end

    return
end

disp(['Reading data from file: ' infile])
dat = load(infile,'out');
dat = dat.out;

if quality  % only use data with rq_flag<=1
    qid = find(dat.rq_flag >= 1 && dat.rq_flag<=3);
    dat.dnum = dat.dnum(qid); 
    dat.deltat = dat.deltat(qid);
    dat.rux = dat.rux(qid); 
    dat.rvy = dat.rvy(qid);
    dat.x = dat.x(qid); 
    dat.y = dat.y(qid);
    dat.stream = dat.stream(qid);
end

% select the data that are within the period of interest
qid = find((dat.dnum >= datenum(startdate)).*(dat.dnum <= datenum(enddate)));
dat.dnum = dat.dnum(qid);
dat.deltat = dat.deltat(qid);
dat.rux = dat.rux(qid);
dat.rvy = dat.rvy(qid);
dat.x = dat.x(qid);
dat.y = dat.y(qid);
dat.stream = dat.stream(qid);

% Initialization
szdat = length(dat.dnum);
streams = unique(dat.stream);
n_streams = length(streams);

out.xy = zeros(2*szdat,2);
out.uv = zeros(2*szdat,2);
out.dnum = zeros(2*szdat,1);
out.deltat = zeros(2*szdat,1);
out.stream = char(ones(2*szdat,1)*'z');

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
    
    [inpath, inname, inext] = fileparts(infile);
    if(filter_noise)
        outfile = fullfile(outdir, [inname(1:end-4) '_deformation_test_min_' num2str(min_size_stream) '_eps_00' num2str(min_def*100) '_level_' num2str(max_level) '.mat']);
    else
        outfile = fullfile(outdir, [inname(1:end-4) '_deformation_test_min_' num2str(min_size_stream) '_no_filter.mat']);
    end
end

% start loop
% Calculate derivatives and strain rates from lagriangian cells
% (triangles).
% Process each stream seperatly.
% Calculation of area A and derivatives after:
% RADARSAT GEOPHYSICAL PROCESSOR SYSTEM, DATA USER???S HANDBOOK (Version 1.0)
% or Kwok et al. (2008), Variability of sea ice simulations assessed with RGPS kinematics, JGR, doi:10.1029/2008JC004783

jstart = 1;
jend=0;
jjstart =1;
first_stream=1;
last_stream=n_streams;
if(test==1)
    first_stream=2;
    last_stream=2;
end
for i=first_stream:last_stream
    
    stream = streams(i);
    id = find(dat.stream == stream);
    
    disp(['Stream ' stream ': ' num2str(i) '/' num2str(n_streams)])
    
    dnum = dat.dnum(id);
    xy = [dat.x(id)'*1000, dat.y(id)'*1000]; % in m
    rux = dat.rux(id).*86400; % in m/day 
    rvy = dat.rvy(id).*86400; % in m/day
    deltat = dat.deltat(id);
    
    min_dnum=min(dnum);
    min_date=min_dnum ;
    if(test)
        min_date=min_date + 65.95 ;%+ 62.9; %+ 83.0;%+ 68.9; %
    end
    % split dataset in orbits (overflights) (0.02days=29min interval)
    dataleft = 1;
    countorb = 0;
    while dataleft
        tmpstart = min(dnum);
        tmpend = tmpstart+0.02;
        tmpid = find(dnum>=tmpstart & dnum<tmpend);
        tdnum = dnum(tmpid);
        dnum(tmpid) = [];
        tdeltat = deltat(tmpid);
        deltat(tmpid) = [];
        txy = xy(tmpid,:);
        xy(tmpid,:) = [];
        trux = rux(tmpid);
        rux(tmpid)  = [];
        trvy = rvy(tmpid);
        rvy(tmpid) = [];
        
        if(tmpstart<min_date)
            continue
        else
            tmpstart-min_dnum;
        end
        
        if isempty(dnum); dataleft = 0; end
        
        % split orbit in sections with ~ the same delta T (-> referring to the same
        % n+1 orbit; different delta T causes artifacts in the triangulation otherwise)
        dataleft2 = 1;
        while dataleft2
            deltstart = min(tdeltat);
            deltend = deltstart+0.02;
            tmpid2 = find(tdeltat>=deltstart & tdeltat<deltend);
            tdeltat2 = mean(tdeltat(tmpid2));
            tdeltat(tmpid2) = [];
            tdnum2 = mean(tdnum(tmpid2));
            tdnum(tmpid2) = [];
            txy2 = txy(tmpid2,:);
            txy(tmpid2,:) = [];
            trux2 = trux(tmpid2);
            trux(tmpid2)  = [];
            trvy2 = trvy(tmpid2);
            trvy(tmpid2) = [];
            
            if isempty(tdeltat); dataleft2 = 0; end
            
            % Limit the domain
            x=txy2(:,1);
            y=txy2(:,2);
            
            
            % Only select te nodes that are within limit_domain
            if(~isempty(limit_domain))
                nb_subdomains=length(limit_domain);
                inside_subdomain=zeros(length(x),1);
                
                x_target=x/1000;
                y_target=y/1000;
                for j=1:nb_subdomains
                    sign_cross_product=zeros(length(x),1);
                    tmp_subdomain=limit_domain(j);
                    nb_edges_polygon=length(tmp_subdomain.x)-1;
                    for k=1:nb_edges_polygon,
                        x1=tmp_subdomain.x(k);
                        x2=tmp_subdomain.x(k+1);
                        y1=tmp_subdomain.y(k);
                        y2=tmp_subdomain.y(k+1);
                        sign_cross_product=sign_cross_product+sign((x_target-x1).*(y2-y1)-(y_target-y1).*(x2-x1));
                    end
                    inside_subdomain=inside_subdomain+(abs(sign_cross_product)==nb_edges_polygon);
                end
                
                ind_selected=find(inside_subdomain);
                
                txy2 =txy2(ind_selected,:);
                trux2=trux2(ind_selected);
                trvy2=trvy2(ind_selected);
            end
            
            
%             if(length(limit_domain)==4)
%                 ind_selected=find((x>=limit_domain(1)).*(x<=limit_domain(2)).*(y>=limit_domain(3)).*(y<=limit_domain(4)));
%                 if(length(ind_selected)==0)
%                     continue;
%                 end
%                 txy2 =txy2(ind_selected,:);
%                 trux2=trux2(ind_selected);
%                 trvy2=trvy2(ind_selected);
%             end
            
            % remove double data points (very important! DelaunyTri gives
            % wrong results with double data points!)
            [txy2, tmpid2] = unique(txy2,'rows');
            trux2 = trux2(tmpid2);
            trvy2 = trvy2(tmpid2);
            
            n2 = length(tmpid2);
            if(compute_deformation)
                if n2>min_size_stream
                    %disp(['n2: ' num2str(n2)])
                    
                    [ntri,ttdnum,ttdeltat,area,xy_tricorner,uv,uv_tricorner,dudx,dudy,dvdx,dvdy,id_stream,quality_index]=def_from_drift(txy2,tdeltat2,tdnum2,trux2,trvy2,stream,min_area,max_area,min_minang,max_long_side,showplot,filter_noise,min_def,max_level);
                    
                    if(ntri>0)
                        % add to output structure
                        jend = jstart+ntri-1;
                        out.dnum(jstart:jend) = ttdnum;
                        out.deltat(jstart:jend) = ttdeltat;
                        out.area(jstart:jend) = area/1e6; % km^2
                        out.xy_tricorner(jstart:jend,:,:) = xy_tricorner/1000; % in km
                        out.xy(jstart:jend,:) = mean(xy_tricorner,2)/1000; % in km
                        out.uv(jstart:jend,:) = uv./1000; % in km/day
                        out.uv_tricorner(jstart:jend,:,:) = uv_tricorner/1000;
                        out.dudx(jstart:jend) = dudx; out.dudy(jstart:jend) = dudy; % in /day
                        out.dvdx(jstart:jend) = dvdx; out.dvdy(jstart:jend) = dvdy; % in /day
                        out.stream(jstart:jend) = id_stream;
                        out.size_datatake(jstart:jend) = ntri;
                        
                        out.quality_index_tot(jjstart)=quality_index;
                        
                        jstart = jend+1;
                        jjstart=jjstart+1;
                    end  
                end % / n observations in orbit
            else
                % add to output structure
                jend = jstart+n2-1;
                out.dnum(jstart:jend) = tdnum2;
                out.deltat(jstart:jend) = tdeltat2;
                out.xy(jstart:jend,:) = txy2./1000; % in km/day
                out.uv(jstart:jend,1) = trux2./1000; % in km/day
                out.uv(jstart:jend,2) = trvy2./1000; % in km/day
                out.stream(jstart:jend) = stream;
                                
                jstart = jend+1;
                jjstart=jjstart+1;
            end
        end % /while different delta T
        countorb = countorb+1;
    end % /orbits
end % /streams


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
