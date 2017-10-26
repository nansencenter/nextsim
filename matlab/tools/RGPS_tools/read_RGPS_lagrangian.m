function out = read_RGPS_lagrangian(year,doy,productl,stream,single_file)
% 
%  Purpose:
%          read_RGPS_lagrangian - opens the lagrangian ice motion product and prints to an
%                   output file the data in the metadata record, the image
%                   description record and ice motion data for the 
%                   first 3 cells in the product
%  
%  Syntax: rd30_lim 
%  Input:  year and day of year 'doy'
%          productl : 'Motion', 'Deformation' [default 'Motion']
%          stream : data stream id, e.g. 'a'
%          single_file: do not read all files for the given doy but only
%                       the 'shortest', 'longest', or matching <fileid>
%  Author: original IDL version rd30_lim: Lisa Nguyen, 03/25/00
%          Matlab conversion: Gunnar Spreen, 2009-02-20
%          more options: Gunnar Spreen, 2012-01-20
%          added summer data, 2012-03
%


fformat = 'ieee-be';

% if different from 0 read only one file
if nargin < 5
    single_file = 'shortest'; % 'shortest':  only read the file with shortest time span of observations
end                           % 'longest': only read the file with longest time span of observations

if nargin < 4
    region_id = 'a' ;  % Data stream ID
else
    region_id = stream ;  % Data stream ID
end

if nargin < 3
    product='LP' ; % lagrangian ice motion
    productl = 'Motion';
else
    switch productl
        case 'Motion'
            product = 'LP';
        case 'Deformation'
            product = 'DP';
            %disp('Only motion implemented at the moment!')
            %return
        otherwise
            disp(['Product ' productl ' not defined!'])
            return
    end
end


% --- end user variables -----
disp(['Using data stream region ' region_id])

if ~ischar(year); year = num2str(year); end
doy = sprintf('%03i',doy);

define_default_data_path;

tmp_dir = [indir 'RGPS_ice_drift/Lagrangian/Winter/' year '/' productl '/'];
file = [tmp_dir,'R1001',region_id,year(3:4),doy,'*.',product] ;

files = dir(file);
nfiles = length(files);

if nfiles == 0; % look for summer data
    tmp_dir = [indir 'RGPS_ice_drift/Lagrangian/Summer/' year '/' productl '/'];
    file = [tmp_dir,'R100*',region_id,year(3:4),doy,'*.',product] ;
    files = dir(file);
    nfiles = length(files);
end

disp([num2str(nfiles) ' files found for date ' year ', ' doy])
if single_file ~= 0
switch single_file
    case 'longest'
        disp('Only reading the file with the longest data record')
        for i=1:nfiles
            bytes(i) = files(i).bytes;
        end
        [~, ind] = max(bytes);
        if ind ~= nfiles
            if bytes(nfiles) ~= bytes(ind)
                disp('Error occured. Check file checking!')
                return
            end
        end
        files = files(nfiles);
        nfiles = 1;
    case 'shortest'
        disp('Only reading the file with the shortest data record')
        for i=1:nfiles
            bytes(i) = files(i).bytes;
        end
        [~, ind] = min(bytes);
        if ind ~= 1
            disp('Error occured. Check file checking!')
            return
        end
        files = files(1);
        nfiles = 1;
    otherwise
        disp(['Only reading file ' single_file])
        for i=1:nfiles
            if strcmp(files(i).name,single_file)
                break
            end
        end
        files = files(i);
        nfiles = 1; 
end
end

for i=1:nfiles
    if nfiles > 1, disp(['Reading file ' num2str(i)]); end
    f = fullfile(tmp_dir,files(i).name);
    disp(f)
    id = fopen(f,'r');
    
    % Lagrangian Ice Motion Product: Metadata Record Contents
    switch productl
        case 'Motion'
        % RGPS product identifier
        out(i).meta.pid = transpose(fread(id,24,'*char')) ;
        % Description of this product
        out(i).meta.prod_desc = transpose(fread(id,40,'*char')) ;
        % Number of images used in the creation of this product
        out(i).meta.n_images = fread(id,1,'short',fformat);
        % Number of trajectories in this product
        out(i).meta.n_trajectories = fread(id,1,'long',fformat) ;
        % Product Type
        out(i).meta.prod_type = transpose(fread(id,8,'*char')) ;
        % Product creation year/time
        out(i).meta.create_year = fread(id,1,'short',fformat) ;
        out(i).meta.create_time = fread(id,1,'double',fformat) ;
        % Season start year/time
        out(i).meta.season_start_year = fread(id,1,'short',fformat);
        out(i).meta.season_start_time = fread(id,1,'double',fformat);
        % Season end year/time
        out(i).meta.season_end_year = fread(id,1,'short',fformat);
        out(i).meta.season_end_time = fread(id,1,'double',fformat);
        % Software version used to create this product
        out(i).meta.sw_version = transpose(fread(id,12,'*char')) ;
        % Northwest Lat/Long of initial datatake
        out(i).meta.n_w_lat = fread(id,1,'float',fformat);
        out(i).meta.n_w_lon = fread(id,1,'float',fformat);
        % ; Northeast Lat/Long of initial datatake
        out(i).meta.n_e_lat = fread(id,1,'float',fformat);
        out(i).meta.n_e_lon = fread(id,1,'float',fformat);
        % Southwest Lat/Long of initial datatake
        out(i).meta.s_w_lat = fread(id,1,'float',fformat);
        out(i).meta.s_w_lon = fread(id,1,'float',fformat);
        % Southeast Lat/Long of initial datatake
        out(i).meta.s_e_lat = fread(id,1,'float',fformat);
        out(i).meta.s_e_lon = fread(id,1,'float',fformat);

        for ii = 1:out(i).meta.n_images
            %   readu,unitr,image_id,image_year,image_time,map_x,map_y
            %fprintf(1,['Image number ',num2str(ii),':']);
            % ASF image identifier
            out(i).meta.image(ii).id = transpose(fread(id,16,'*char')); 
            % Image center year/time
            out(i).meta.image(ii).year = fread(id,1,'short',fformat) ;
            out(i).meta.image(ii).time = fread(id,1,'double',fformat) ;
            % Image center location
            out(i).meta.image(ii).map_x = fread(id,1,'double',fformat);
            out(i).meta.image(ii).map_y = fread(id,1,'double',fformat);   
        end
        
        % =======================================================
        %  GRIDPOINT/TRAJECTORY DESCRIPTION DATA
        % =======================================================
        % Read the Lagrangian Ice Motion Product: Gridpoint/Trajectory Description Data
        %
        % The gridpoint description data contain one record for each gridpoint
        % and its trajectory.
        %
        % printf,unitw,""
        % printf,unitw,"LAGRANGIAN ICE MOTION PRODUCT: Ice Motion Data"
        % printf,unitw,"----------------------------------------------"

        %ii=0;
        %while ~feof(id)
        for ii=1:out(i).meta.n_trajectories
            %ii = ii+1;
            % Grid point identifier
            out(i).trajectories(ii).gpid = fread(id,1,'long',fformat);
            % Birth and Death year/time of gridpoint
            out(i).trajectories(ii).birth_year = fread(id,1,'short',fformat);
            out(i).trajectories(ii).birth_time = fread(id,1,'double',fformat);
            out(i).trajectories(ii).death_year = fread(id,1,'short',fformat);
            out(i).trajectories(ii).death_time = fread(id,1,'double',fformat);
            % Number of observations in trajectory
            out(i).trajectories(ii).n_obs = fread(id,1,'long',fformat);

            % printf,unitw,"  GRIDID   N_OBS   B_YEAR    B_TIME   D_YEAR    D_TIME"
            % printf,unitw,gpid,n_obs,birth_year,birth_time,death_year,death_time,$
            %        format="(i7,i8,i9,f12.5,i7,f12.5)"
            %
            % Read the each information of each observation in trajectory
            %
            for iii=1:out(i).trajectories(ii).n_obs
                % Year/Time of observation
                out(i).trajectories(ii).year(iii) = fread(id,1,'short',fformat);
                out(i).trajectories(ii).time(iii) = fread(id,1,'double',fformat);
                % Map location of observation
                out(i).trajectories(ii).x_map(iii) = fread(id,1,'double',fformat);
                out(i).trajectories(ii).y_map(iii) = fread(id,1,'double',fformat);
                % Quality Flag of observation
                out(i).trajectories(ii).q_flag(iii) = fread(id,1,'short',fformat);
            end
            [lon,lat]=xyssmi1p((out(i).trajectories(ii).x_map(:))*1000,(out(i).trajectories(ii).y_map(:))*1000,-45);
            m_proj('Stereographic','lon',-45,'lat',90,'radius',60);
            [x,y]=m_ll2xy(lon,lat);
            x=(x*6378.273);
            y=(y*6378.273);
            out(i).trajectories(ii).x_map(:) = x;
            out(i).trajectories(ii).y_map(:) = y;
            %disp(100*i/nfiles)
        end        
        
        
        case 'Deformation'
            % RGPS product identifier
            out(i).meta.idf_id = transpose(fread(id,24,'*char')) ;
            % Description of this product
            out(i).meta.prod_desc = transpose(fread(id,40,'*char')) ;
            % Number of cells in this product
            out(i).meta.n_cells = fread(id,1,'long',fformat);
            % Product creation year/time
            out(i).meta.create_year = fread(id,1,'short',fformat) ;
            out(i).meta.create_time = fread(id,1,'double',fformat) ;
            % Season start year/time
            out(i).meta.season_start_year = fread(id,1,'short',fformat);
            out(i).meta.season_start_time = fread(id,1,'double',fformat);
            % Season end year/time
            out(i).meta.season_end_year = fread(id,1,'short',fformat);
            out(i).meta.season_end_time = fread(id,1,'double',fformat);
            % Software version used to create this product
            out(i).meta.sw_version = transpose(fread(id,12,'*char')) ;
            % Northeast Lat/Long of initial datatake
            out(i).meta.n_e_lat = fread(id,1,'float',fformat);
            out(i).meta.n_e_lon = fread(id,1,'float',fformat);
            % ; Northwest Lat/Long of initial datatake
            out(i).meta.n_w_lat = fread(id,1,'float',fformat);
            out(i).meta.n_w_lon = fread(id,1,'float',fformat);
            % Southwest Lat/Long of initial datatake
            out(i).meta.s_w_lat = fread(id,1,'float',fformat);
            out(i).meta.s_w_lon = fread(id,1,'float',fformat);
            % Southeast Lat/Long of initial datatake
            out(i).meta.s_e_lat = fread(id,1,'float',fformat);
            out(i).meta.s_e_lon = fread(id,1,'float',fformat);
            
            for ii=1:out(i).meta.n_cells
                %ii = ii+1;
                % Cell identifier
                out(i).cells(ii).cellid = fread(id,1,'long',fformat);
                % Birth and Death year/time of cell
                out(i).cells(ii).birth_year = fread(id,1,'short',fformat);
                out(i).cells(ii).birth_time = fread(id,1,'double',fformat);
                % Number of observations
                out(i).cells(ii).n_obs = fread(id,1,'short',fformat);
                
                % printf,unitw,"  GRIDID   N_OBS   B_YEAR    B_TIME   D_YEAR    D_TIME"
                % printf,unitw,gpid,n_obs,birth_year,birth_time,death_year,death_time,$
                %        format="(i7,i8,i9,f12.5,i7,f12.5)"
                %
                % Read the each information of each observation in trajectory
                %
                for iii=1:out(i).cells(ii).n_obs
                    % Year/Time of observation
                    out(i).cells(ii).year(iii) = fread(id,1,'short',fformat);
                    out(i).cells(ii).time(iii) = fread(id,1,'double',fformat);
                    % Map location of observation
                    out(i).cells(ii).x_map(iii) = fread(id,1,'double',fformat);
                    out(i).cells(ii).y_map(iii) = fread(id,1,'double',fformat);
                    % Displacements
                    out(i).cells(ii).x_disp(iii) = fread(id,1,'double',fformat);
                    out(i).cells(ii).y_disp(iii) = fread(id,1,'double',fformat);
                    % Areas
                    out(i).cells(ii).c_area(iii) = fread(id,1,'float',fformat);
                    out(i).cells(ii).d_area(iii) = fread(id,1,'float',fformat);
                    % Time interval
                    out(i).cells(ii).dtp(iii) = fread(id,1,'float',fformat);
                    % Derivatives of the displacement
                    out(i).cells(ii).dudx(iii) = fread(id,1,'float',fformat);
                    out(i).cells(ii).dudy(iii) = fread(id,1,'float',fformat);
                    out(i).cells(ii).dvdx(iii) = fread(id,1,'float',fformat);
                    out(i).cells(ii).dvdy(iii) = fread(id,1,'float',fformat);
                end
                [lon,lat]=xyssmi1p((out(i).cells(ii).x_map(:))*1000,(out(i).cells(ii).y_map(:))*1000,-45);
                m_proj('Stereographic','lon',-45,'lat',90,'radius',60);
                [x,y]=m_ll2xy(lon,lat);x=(x*6378.273);y=(y*6378.273);
                out(i).cells(ii).x_map(:) = x;
                out(i).cells(ii).y_map(:) = y;
            end
        otherwise
            disp(['Product ' productl ' not defined!'])
            return
    end
    
    fclose(id);
end % for files loop