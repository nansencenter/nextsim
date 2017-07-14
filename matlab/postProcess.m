function postProcess(directory,skip,outputfile)

% function postProcess(directory,skip,outputfile)
%
% directory:    The directory where the output is stored
% skip:         Take only every n-th output
% outputfile:   Save output to mat file to use with volumePlot.m

global mppfile
mppfile = which('NpsNextsim.mpp');
%addpath('~/src/nextsim/matlab/')

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
% R = 6378.273;
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

% [x_ICE, y_ICE] = m_ll2xy(maskll(:,1),maskll(:,2));
[x_ICE, y_ICE] = mapx_forward(mppfile, maskll(:,1)',maskll(:,2)');

% Now points south of the DRA north of Greenland, Canada, and in the
% Beaufort Sea
maskll = [ -141.0000   70.0000
           -141.0000   80.0000
           -130.0000   80.0000
            -60.0000   86.5800
            -15.0000   87.0000
              8.0000   84.0800
            -12.2516   81.5707
            -47.9017   80.9993
            -73.7776   83.7842
           -103.2309   79.6701
           -121.9610   77.1352
           -128.5021   68.9701
           -140.9908   68.9919 ];
       
[x_GCB, y_GCB] = mapx_forward(mppfile, maskll(:,1)',maskll(:,2)');

% And points south of the DRA, north of Greenland and Canada only
maskll = [ -130.0000   80.0000
            -60.0000   86.5800
            -15.0000   87.0000
              8.0000   84.0800
            -12.2516   81.5707
            -47.9017   80.9993
            -73.7776   83.7842
           -103.2309   79.6701
           -121.9610   77.1352 ];
       
[x_GC, y_GC] = mapx_forward(mppfile, maskll(:,1)',maskll(:,2)');

% Then we load all the data
files = dir([directory '/field_*.bin']);

% List and sort the available output steps
steps = zeros(length(files),1);
for i=1:length(files)
    steps(i) = str2double(files(i).name(7:end-4));
end
steps = sort(steps);
steps(isnan(steps)) = []; % Throwing out the _init file
if steps(end) > steps(end-1)+1
    steps(end) = []; % Throwing out the 1000 file - if it doesn't belong
end

% Initialize for speed
t = 0.*(steps(1):skip:steps(end));

vice     = t;
vice_DRA = t;
vice_ICE = t;

area       = t;
area_SSMI  = t;
area_AMSRE = t;

extent       = t;
extent_SSMI  = t;
extent_AMSRE = t;

vsnow = t;

k=0;
for i = steps(1):skip:steps(end)
    k=k+1;
    
    [simul_out, elementx, elementy, element_lat] = bin2simul_out(directory, i);
    
    mask_DRA = inpolygon(elementx,elementy,x_DRA,y_DRA);
    mask_ICE = inpolygon(elementx,elementy,x_ICE,y_ICE);
    mask_GCB = inpolygon(elementx,elementy,x_GCB,y_GCB);
    mask_GC = inpolygon(elementx,elementy,x_GC,y_GC);
    mask_SSMI  = element_lat <= 87.24;
    mask_AMSRE = element_lat <= 89.24;

    t(k) = simul_out.current_time;

    vice(k)   = simul_out.h'*simul_out.surface;
    vice_DRA(k)   = simul_out.h(mask_DRA)'*simul_out.surface(mask_DRA);
    vice_ICE(k)   = simul_out.h(mask_ICE)'*simul_out.surface(mask_ICE);
    vice_GCB(k)   = simul_out.h(mask_GCB)'*simul_out.surface(mask_GCB);
    vice_GC(k)   = simul_out.h(mask_GC)'*simul_out.surface(mask_GC);

    area(k)       = simul_out.c'*simul_out.surface;
    area_SSMI(k)  = simul_out.c(mask_SSMI)'*simul_out.surface(mask_SSMI);
    area_AMSRE(k) = simul_out.c(mask_AMSRE)'*simul_out.surface(mask_AMSRE);
    
    extent(k)       = (simul_out.c >= 0.15)'*simul_out.surface;
    extent_SSMI(k)  = (simul_out.c(mask_SSMI) >= 0.15)'*simul_out.surface(mask_SSMI);
    extent_AMSRE(k) = (simul_out.c(mask_AMSRE) >= 0.15)'*simul_out.surface(mask_AMSRE);

    vsnow(k)     = simul_out.hs'*simul_out.surface;
end

save(outputfile, 'directory', 'vice', 'area', 'extent', 'vsnow', 't', ...
    'vice_ICE', 'vice_DRA', 'vice_GCB', 'vice_GC', ...
    'area_SSMI', 'area_AMSRE', 'extent_SSMI', 'extent_AMSRE')

end

function [simul_out, elementx, elementy, element_lat] = bin2simul_out(dir,step)

% Load the data
[mesh_out, data_out] = neXtSIM_bin_revert(dir, '', step);

% Calculate coordinates for elements
% reshape
var_mx=mesh_out.Nodes_x(mesh_out.Elements);
var_my=mesh_out.Nodes_y(mesh_out.Elements);

[nr,nc]= size(var_mx);
Ne=nr/3;
Nn=length(mesh_out.Nodes_x);
elementx=mean(reshape(var_mx,[3,Ne]));
elementy=mean(reshape(var_my,[3,Ne]));

% R = 6378.273;
% [~,element_lat] = m_xy2ll(elementx/R, elementy/R);
% mppfile = '../data/NpsNextsim.mpp';
global mppfile;
[~, element_lat] = mapx_inverse(mppfile, elementx, elementy);

% Put things in the simul_out struct
simul_out.current_time  = data_out.Time;
simul_out.h             = data_out.Thickness;
simul_out.c             = data_out.Concentration;
simul_out.surface       = data_out.Element_area;
simul_out.hs            = data_out.Snow;
try % and see if we're using thin ice
    simul_out.h    = simul_out.h  + data_out.Thin_ice;
    simul_out.hs   = simul_out.hs + data_out.Snow_thin_ice;
    simul_out.c    = simul_out.c  + data_out.Concentration_thin_ice;
end

end
