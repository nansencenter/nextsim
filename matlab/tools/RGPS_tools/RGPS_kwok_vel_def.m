function RGPS_kwok_vel_def(year,productl)
% Produces RGPS lagrangian ice motion and deformation dataset.
%
% v1.0 Jan 2012, Gunnar Spreen

%clear all
% startdate:
%year = 1997;
summer = 0; % use summer and not winter data

rgps_file_length = 'longest'; % only read the 'shortest', 'longest', or 
                               % <fileid> RGPS file for a given year and doy

if(strcmp(productl,'Trajectories')||strcmp(productl,'Motion'))
    product_Lagarngian='Motion';
    product='LP' ; % lagrangian ice motion
elseif(strcmp(productl,'Deformation'))
    product='DP' ; % lagrangian ice motion
else
    error('Do you want the Deformation or the Motion field??')
end
   
define_default_data_path;

remove_problematic_stream_a=0; % stream a could be problematic (Ask Gunnar Spreen !!!!)

streamremove = {};
switch year
    case 1996
        doy = 312;
    case 1997
        doy = 306;
        if summer; doy = 138; end
    case 1998
        doy = 301;
        if summer; doy = 130; end
    case 1999
        doy = 305;
        if summer; doy = 128; end
    case 2000
        doy = 309;
        if(strcmp(product,'LP') && remove_problematic_stream_a)
            streamremove = {'a'}; % problem with stream a
        end
    case 2001
        if summer; doy = 135;
        else
            doy = 309;
            if(strcmp(product,'LP') && remove_problematic_stream_a)
                streamremove = {'a'}; % problem with stream a
            end
        end
    case 2002
        doy = 361;
        if(strcmp(product,'LP') && remove_problematic_stream_a)
            streamremove = {'a'}; % problem with stream a
        end
        if summer; doy = 136; end
    case 2003
        doy = 338;
        if(strcmp(product,'LP') && remove_problematic_stream_a)
            streamremove = {'a'}; % problem with stream a
        end
    case 2004
        doy = 315;
        if summer; doy = 132; end
        if(strcmp(product,'LP') && remove_problematic_stream_a)
            streamremove = {'a'}; % problem with stream a
        end
    case 2005
        doy = 333;
        if summer; doy = 135; end
        if(strcmp(product,'LP') && remove_problematic_stream_a)
            streamremove = {'a'}; % problem with stream a
        end
    case 2006
        doy = 337;
        if summer; doy = 139; end
        if(strcmp(product,'LP') && remove_problematic_stream_a)
            streamremove = {'a'}; % problem with stream a
        end
    case 2007
        doy = 335;
        if summer; doy = 134; end
        if(strcmp(product,'LP') && remove_problematic_stream_a)
            streamremove = {'a'}; % problem with stream a
        end
end

% ----------------------------

% -------------------
% read RGPS to memory
% -------------------
if ~exist('doy','var')
    doy = date2doy(datenum(year,month,day));
    doy = doy(2);
else
    [day,month] = doy2date(doy,year)
end
syear = num2str(year);
sdoy = num2str(doy);

if summer
    indir = [indir 'RGPS_ice_drift/Lagrangian/Summer/' syear '/Motion/'];
else
    indir = [indir 'RGPS_ice_drift/Lagrangian/Winter/' syear '/Motion/'];
end
files = dir([indir,'R100*' syear(3:4) '*' sdoy,'*.LP']);
streams = cell(1,length(files));
for i = 1:length(files) % count number of streams
    streams{i} = files(i).name(6);
end
streams = unique(streams);
if ~isempty(streamremove)
    keepid = ~strcmp(streams, streamremove);
    streams = streams(keepid);
end
n_streams = length(streams);
startyears = zeros(n_streams,1);
starttimes =  zeros(n_streams,1);
endyears = zeros(n_streams,1);
endtimes =  zeros(n_streams,1);
usedstreams = cell(1,n_streams);
stream_info = cell(1,n_streams);

disp(['Reading RGPS data for ' num2str(n_streams) ' streams'])
j = 0;
for i = 1:n_streams
    stream = streams{i};
    tmp = read_RGPS_lagrangian(year,doy,product_Lagarngian,stream,rgps_file_length);
    
    switch productl
        case 'Trajectories'
            n_objects = tmp.meta.n_trajectories;
        case 'Motion'
            n_objects = tmp.meta.n_trajectories;
        case 'Deformation'
            n_objects = tmp.meta.n_cells;
        otherwise
            disp(['Product ' productl ' not defined!'])
            return
    end
    
    if n_objects ~= 0
        j = j+1;
        usedstreams{j} = stream;
        tmp.name = stream;
        stream_info{j} = tmp;
        
        lat = [stream_info{j}.meta.n_w_lat, stream_info{j}.meta.n_e_lat, stream_info{j}.meta.s_e_lat, ...
            stream_info{j}.meta.s_w_lat, stream_info{j}.meta.n_w_lat];
        lon = [stream_info{j}.meta.n_w_lon, stream_info{j}.meta.n_e_lon, stream_info{j}.meta.s_e_lon, ...
            stream_info{j}.meta.s_w_lon, stream_info{j}.meta.n_w_lon];
        startyears(j) = tmp.meta.season_start_year;
        starttimes(j) = tmp.meta.season_start_time;
        endyears(j) = tmp.meta.season_end_year;
        endtimes(j) = tmp.meta.season_end_time;

        % show map of data regions
        if i==1
            figure(1)
            worldmap([65 90], [180 360]);
            geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
        end
        geoshow(lat,lon)
    end
end
% shorten records if some streams were empty
n_streams = j;
streams = usedstreams(1:j);
stream_info=stream_info(1:j);
startyears = startyears(1:j);
starttimes = starttimes(1:j);
endyears = endyears(1:j);
endtimes = endtimes(1:j);

startyear = min(startyears);
starttime = min(starttimes);
[startday,startmonth] = doy2date(starttime,startyear);
startdnum = datenum(startyear,startmonth,startday);
startdate = datestr(startdnum,29);
endyear = max(endyears);
endtime = max(endtimes);
[endday,endmonth] = doy2date(endtime,endyear);
enddnum = datenum(endyear,endmonth,endday);
enddate = datestr(enddnum,29);
ndays = ceil(enddnum - startdnum + 2);

switch productl
    case 'Trajectories'
        out = stream_info;
    case 'Motion'
        out = struct('dnum',zeros(1,ndays*13000), 'rq_flag',zeros(1,ndays*13000), ...
            'deltat',zeros(1,ndays*13000), 'x',zeros(1,ndays*13000), 'y',zeros(1,ndays*13000), ...
            'rux',zeros(1,ndays*13000), 'rvy',zeros(1,ndays*13000), ...
            'stream',char(zeros(1,ndays*13000)) );
        n_objects = tmp.meta.n_trajectories;
    case 'Deformation'
        out = struct('dnum',zeros(1,ndays*13000), ...
            'deltat',zeros(1,ndays*13000), 'x',zeros(1,ndays*13000), 'y',zeros(1,ndays*13000), ...
            'opening_rate',zeros(1,ndays*13000), 'area',zeros(1,ndays*13000),'init_area',zeros(1,ndays*13000),...
            'dudx',zeros(1,ndays*13000), 'dudy',zeros(1,ndays*13000), ...
            'dvdx',zeros(1,ndays*13000), 'dvdy',zeros(1,ndays*13000), ...
            'stream',char(zeros(1,ndays*13000)) );      
        n_objects = tmp.meta.n_cells;
    otherwise
        disp(['Product ' productl ' not defined!'])
        return
end
% allocate memory for output structure

j=0;
for i = 1:n_streams
    disp([num2str(i) '/' num2str(n_streams) ' streams'])
    
    switch productl
        case 'Trajectories'
            continue
        case 'Motion'
            n_objects = stream_info{i}.meta.n_trajectories;
            type_objects = 'trajectories';
        case 'Deformation'
            n_objects = stream_info{i}.meta.n_cells;
            type_objects = 'cells';
        otherwise
            disp(['Product ' productl ' not defined!'])
            return
    end
    for ii = 1:n_objects
        if mod(ii,500) == 0
            if exist('tStart','var'); toc(tStart); end
            disp([num2str(ii) '/' num2str(n_objects) ...
                ' ' type_objects ', ' num2str(i) '/' num2str(n_streams) ' streams'])
            tStart=tic;
        end
        
        switch productl
        case 'Motion'
            tmp_obj = stream_info{i}.trajectories(ii);
            first_obs=2;
        case 'Deformation'
            tmp_obj = stream_info{i}.cells(ii);
            first_obs=1;
            if(tmp_obj.n_obs>=1) init_area = tmp_obj.c_area(1); end
        otherwise
            disp(['Product ' productl ' not defined!'])
            return
        end
        
        skipcount = 0;
        for iii = first_obs:tmp_obj.n_obs
            tmp = tmp_obj;
            date = datenum(tmp.year(iii  ), 0, tmp.time(iii));
            date_before = date;
            process = 1;
            if(first_obs>1)
                date_before = datenum(tmp.year(iii-1), 0, tmp.time(iii-1)); 
                if isempty(tmp.year(iii-1)) | tmp.year(iii-1) < 1000 %#ok<OR2>
                    skipcount = skipcount + 1;
                    process = 0;
                end
            end
            if process
                j = j+1;

                if date_before < startdnum
                    disp('wrong date')
                    return
                end
                switch productl
                    case 'Motion'
                        out.rq_flag(j) = tmp.q_flag(iii-1);
                        out.deltat(j) = date - date_before;
                        % RGPS velocities
                        out.dnum(j) = date;
                        out.rux(j) = (tmp.x_map(iii) - tmp.x_map(iii-1))*1000/out.deltat(j)/24/60/60;
                        out.rvy(j) = (tmp.y_map(iii) - tmp.y_map(iii-1))*1000/out.deltat(j)/24/60/60;
                        out.x(j) = tmp.x_map(iii); out.y(j) = tmp.y_map(iii);
                        out.stream(j) = streams{i};
                    case 'Deformation'
                        if(tmp.dtp(iii)==0)
                            % RGPS deformation rates
                            j=j-1;
                        else
                            out.dnum(j) = date;
                            out.deltat(j) = tmp.dtp(iii);
                            out.area(j) = tmp.c_area(iii);
                            out.init_area(j) = init_area;
                            % RGPS deformation rates (/day)
                            out.dudx(j) = tmp.dudx(iii)/out.deltat(j);
                            out.dudy(j) = tmp.dudy(iii)/out.deltat(j);
                            out.dvdx(j) = tmp.dvdx(iii)/out.deltat(j);
                            out.dvdy(j) = tmp.dvdy(iii)/out.deltat(j);
                            out.opening_rate(j) = tmp.d_area(iii)/out.deltat(j);
                            out.x(j) = tmp.x_map(iii); out.y(j) = tmp.y_map(iii);
                            out.stream(j) = streams{i};
                        end
                     otherwise
                        disp(['Product ' productl ' not defined!'])
                        return    
                end
            end % /process observation
        end % /iii obs
        if skipcount ~= 0
            disp(['Faulty record(s) found! Skipped ' num2str(skipcount) ...
                '/' num2str(tmp_obj.n_obs-1) ' obs in stream ' num2str(i) ...
                ', ' type_objects ' ' num2str(ii)])
        end
    end % /ii traject
end % /i streams

switch productl
    case 'Trajectories'
        % cut output array to used length
        outfile = ['RGPS_' startdate '_' enddate '_traj.mat'];
    case 'Motion'
        % cut output array to used length
        out.dnum = out.dnum(1:j); out.rq_flag = out.rq_flag(1:j); out.deltat = out.deltat(1:j);
        out.x = out.x(1:j); out.y = out.y(1:j); out.rux = out.rux(1:j); out.rvy = out.rvy(1:j);
        out.stream = out.stream(1:j);
        outfile = ['RGPS_' startdate '_' enddate '_vel.mat'];
    case 'Deformation'
        % cut output array to used length
        out.dnum = out.dnum(1:j); out.deltat = out.deltat(1:j);
        out.opening_rate = out.opening_rate(1:j); out.area = out.area(1:j);
        out.init_area = out.init_area(1:j);
        out.x = out.x(1:j); out.y = out.y(1:j); 
        out.dudx = out.dudx(1:j); out.dudy = out.dudy(1:j);
        out.dvdx = out.dvdx(1:j); out.dvdy = out.dvdy(1:j);
        out.stream = out.stream(1:j);
        outfile = ['RGPS_' startdate '_' enddate '_def.mat'];
    otherwise
        disp(['Product ' productl ' not defined!'])
        return
end
disp(['Save data to ' outfile])
save(outfile,'out');

beep
disp('End')


