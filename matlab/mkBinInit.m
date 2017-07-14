function mkBinInit(data, lon, lat, fname, meshdir, mppfile)

% A simple script to take any field with associated longitude and latitude
% and turn into initial conditions to be read using --setup.ice-type=binary
%
% Usage:
%
% >> mkBinInit(data, lon, lat, file_name, [mesh_dir, [mppfile]])
%
% data: The data field to be saved (concentration, thickness, or snow thickness)
% lon, lat: The associated geographical coordinates
% filename: Name of the output file (initConc.dat, initThick.dat, or initSnow.dat).
%       To be saved in $NEXTSIMDIR/data/
% meshdir: Directory where the mesh_0 and field_0 files are found
%       Optional, defaults to the current directory
% mppfile: Location of the projection file
%       Optional, defaults to ../data/NpsNextsim.mpp
%
% NB. 
% * In order to have mesh_0 and field_0 files the model must be run for
%       at least one time step using some (any) other initial conditions.
% * When run with --setup.ice-type=binary the model expects to find files
%       named initConc.dat (concentration), initThick.dat (ice thickness), 
%       and initSnow.dat (snow thickness) in $NEXTSIMDIR/data/

% Inputs
if nargin < 6
    mppfile = which('NpsNextsim.mpp');
    if nargin < 5
        meshdir = '';
        if nargin < 4
            error('Not enough input arguments')
        end
    end
end

% Mask out NaNs in input data
mask = ~isnan(data);

% Project
[xdata, ydata] = mapx_forward(mppfile, lon(mask)', lat(mask)');

% Load the mesh and calculate x and y for the elements
mesh_out = neXtSIM_bin_revert(meshdir,'',0);
var_mx = mesh_out.Nodes_x(mesh_out.Elements);
var_my = mesh_out.Nodes_y(mesh_out.Elements);
nr = size(var_mx,1);
Ne = nr/3;
x = mean(reshape(var_mx,[3,Ne]));
y = mean(reshape(var_my,[3,Ne]));

% Interpolate the input data
idata = griddata(xdata, ydata, data(mask), x, y);

% Write to file
fid=fopen(fname,'w');
fwrite(fid,idata,'double');
fclose(fid);

