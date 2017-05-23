% function volICEsat = ICEsat_prep(tICEsat, meshdir)

% The directory where outputfiles from the model are kept. I interpolate
% everything onto the mesh in the first output file.
meshdir = '';

addpath('~/sim/data/CS2_SMOS_v11/');
addpath('~/src/nextsim/matlab');
files = dir('~/sim/data/CS2_SMOS_v11/cs2_smos_ice_thickness_*.nc');
%files = dir('~/sim/data/CS2_SMOS_thickness/cs2_smos_ice_thickness_2010*.nc');

%addpath('/net/sverdrup-1/vol/sim/data/SIT_data/icesat_filled_10prods')
mppfile = which('NpsNextsim.mpp');

% Start with the masks
% First we define the DRA mask
maskll  = [ -15.00 87.00
            -60.00 86.58
            -130.00 80.00
            -141.00 80.00
            -141.00 70.00
            -155.00 72.00
             175.00 75.50
             172.00 78.50
             163.00 80.50
             126.00 78.50
             110.00 84.33
              80.00 84.42
              57.00 85.17
              33.00 83.83
               8.00 84.08
             -15.00 87.00 ];
% m_proj('Stereographic','lon',-45,'lat',90,'radius',60);
% [x_DRA, y_DRA] = m_ll2xy(maskll(:,1),maskll(:,2));
[x_DRA, y_DRA] = mapx_forward(mppfile, maskll(:,1)',maskll(:,2)');

% Then the ICEsat mask
maskll  = [ -47.9017   80.9993
             17.4651   78.9698
             54.8947   80.0657
             92.7812   80.6256
            116.7218   69.8991
           -188.7011   66.4105
           -166.2534   65.6531
           -128.5021   68.9701
           -121.9610   77.1352
           -103.2309   79.6701
            -73.7776   83.7842
            -47.9017   80.9993 ];

[x_ICE, y_ICE] = mapx_forward(mppfile, maskll(:,1)',maskll(:,2)');
% [x_ICE, y_ICE] = m_ll2xy(maskll(:,1),maskll(:,2));
% R = 6378.273;

% Interpolate onto our grid and then calculate the volume (to get a round a
% slight offset due to different grids and masks).
% load topazreducedsplit2
[mesh, data] = neXtSIM_bin_revert(meshdir,'',0);

% reshape
var_mx=mesh.Nodes_x(mesh.Elements);
var_my=mesh.Nodes_y(mesh.Elements);

[nr,nc]= size(var_mx);
Ne=nr/3;
Nn=length(mesh.Nodes_x);
element.x=mean(reshape(var_mx,[3,Ne]));
element.y=mean(reshape(var_my,[3,Ne]));

element.surf = data.Element_area;

mask_ICE = inpolygon(element.x,element.y,x_ICE,y_ICE);
mask_DRA = inpolygon(element.x,element.y,x_DRA,y_DRA);

ncfile = files(1).name;
lon = ncread(ncfile, 'longitude');
lat = ncread(ncfile, 'latitude');

[x, y] = mapx_forward(mppfile, lon(:)', lat(:)');

tCS2_SMOS   = nan(1,length(files));
volCS2_SMOS    = nan(1,length(files));
volminCS2_SMOS = nan(1,length(files));
volmaxCS2_SMOS = nan(1,length(files));
volCS2_SMOS_ICE    = nan(1,length(files));
volminCS2_SMOS_ICE = nan(1,length(files));
volmaxCS2_SMOS_ICE = nan(1,length(files));
volCS2_SMOS_DRA    = nan(1,length(files));
volminCS2_SMOS_DRA = nan(1,length(files));
volmaxCS2_SMOS_DRA = nan(1,length(files));
for k=1:length(files)

    ncfile = files(k).name;

    conc  = ncread(ncfile, 'ice_concentration');
    thick = ncread(ncfile, 'analysis_thickness');
    err   = ncread(ncfile, 'analysis_thickness_err');
    thick(isnan(thick)) = 0;
    err(isnan(thick))   = 0;
    mask = conc>0; %~isnan(conc);
    
    F = scatteredInterpolant(x(mask), y(mask), double(thick(mask))); %.*conc(mask)/100));
    thick_int = F(element.x, element.y);
    
    F.Values = double(err(mask).*conc(mask)/100);
    err_int  = F(element.x, element.y);

    t0 = datenum(str2double(ncfile(24:27)), str2double(ncfile(28:29)), str2double(ncfile(30:31)));
    t1 = datenum(str2double(ncfile(33:36)), str2double(ncfile(37:38)), str2double(ncfile(39:40)));
    tCS2_SMOS(k) = 0.5*(t0+t1);

    volCS2_SMOS(k)     = nansum(thick_int'.*element.surf);
    volminCS2_SMOS(k)  = nansum(max(0, (thick_int-err_int)'.*element.surf));
    volmaxCS2_SMOS(k)  = nansum(max(0, (thick_int+err_int)'.*element.surf));
    
    volCS2_SMOS_ICE(k)    = nansum(thick_int(mask_ICE)'.*element.surf(mask_ICE));
    volminCS2_SMOS_ICE(k) = nansum(max(0, (thick_int(mask_ICE)-err_int(mask_ICE))'.*element.surf(mask_ICE)));
    volmaxCS2_SMOS_ICE(k) = nansum(max(0, (thick_int(mask_ICE)+err_int(mask_ICE))'.*element.surf(mask_ICE)));
    
    volCS2_SMOS_DRA(k)    = nansum(thick_int(mask_DRA)'.*element.surf(mask_DRA));
    volminCS2_SMOS_DRA(k) = nansum(max(0, (thick_int(mask_DRA)-err_int(mask_DRA))'.*element.surf(mask_DRA)));
    volmaxCS2_SMOS_DRA(k) = nansum(max(0, (thick_int(mask_DRA)+err_int(mask_DRA))'.*element.surf(mask_DRA)));
end

save('CS2_SMOS.mat', 'tCS2_SMOS', 'volCS2_SMOS', 'volminCS2_SMOS', 'volmaxCS2_SMOS', ...
    'volCS2_SMOS_ICE', 'volminCS2_SMOS_ICE', 'volmaxCS2_SMOS_ICE', ...
    'volCS2_SMOS_DRA', 'volminCS2_SMOS_DRA', 'volmaxCS2_SMOS_DRA');
