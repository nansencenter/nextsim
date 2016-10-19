%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simul_out_file = '~/experiments/ref_run_750kPa_28days/simul_out_topazreducedsplit2_750kPa_28days_step0.mat'
simul_in_file  = '~/experiments/ref_run_750kPa_28days/simul_in_topazreducedsplit2_750kPa_28days.mat'
restart_dir    = '~/src/nextsim/restart/'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output precission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int    = 'int32';
double = 'float64';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the simul out and simul in file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(simul_in_file)
load(simul_out_file)
mesh = importbamg(simul_out.bamg.mesh, simul_out.bamg.geom);

indx = regexp(simul_out_file, 'step[0-9]+\.mat$');
step_str = simul_out_file(indx+4:end-4);
step = str2double(step_str);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the mesh files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meshf_dat = [restart_dir '/mesh_' step_str '.dat']
meshf_bin = [restart_dir '/mesh_' step_str '.bin']
fid_dat = fopen(meshf_dat, 'w', 'n', 'UTF-8');
fid_bin = fopen(meshf_bin, 'w');

% Connectivity:
% Use [1 3 2] to get the right the orientation
fprintf(fid_dat, 'Elements int\n');
fwrite(fid_bin, numel(mesh.element.num_node), int);
fwrite(fid_bin, mesh.element.num_node(:,[1 3 2])', int);

% Coordinates:
% Convert to lat/lon using m_xy2ll and then back to xy using mapx_forward
mppfile = '../data/NpsNextsim.mpp';
m_proj('Stereographic','lon',-45,'lat',90,'radius',60);
R = 6378.273;
[lon, lat] = m_xy2ll(mesh.node.x/R, mesh.node.y/R);
[x, y] = mapx_forward(mppfile, lon, lat);
fprintf(fid_dat, 'Nodes_x double\n');
fwrite(fid_bin, numel(x), int);
fwrite(fid_bin, x, double); % 

fprintf(fid_dat, 'Nodes_y double\n');
fwrite(fid_bin, numel(y), int);
fwrite(fid_bin, y, double);

fprintf(fid_dat, 'id int');
fwrite(fid_bin, numel(simul_out.node_id), int);
fwrite(fid_bin, simul_out.node_id, int);

fclose(fid_dat);
fclose(fid_bin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now for the field files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fieldf_dat = [restart_dir '/field_' step_str '.dat']
fieldf_bin = [restart_dir '/field_' step_str '.bin']
fid_dat = fopen(fieldf_dat, 'w', 'n', 'UTF-8');
fid_bin = fopen(fieldf_bin, 'w');

fprintf(fid_dat, 'Misc_int int\n');
fwrite(fid_bin, 2, int);
fwrite(fid_bin, [simul_out.p simul_in.flag_fix], int);

fprintf(fid_dat, 'M_dirichlet_flags int\n');
fwrite(fid_bin, numel(simul_in.ind_node_fix_bnd), int);
fwrite(fid_bin, simul_in.ind_node_fix_bnd-1, int); % -1 to start at 0

fprintf(fid_dat, 'M_conc double\n');
fwrite(fid_bin, numel(simul_out.c), int);
fwrite(fid_bin, simul_out.c, double);

fprintf(fid_dat, 'M_thick double\n');
fwrite(fid_bin, numel(simul_out.h), int);
fwrite(fid_bin, simul_out.h, double);

fprintf(fid_dat, 'M_snow_thick double\n');
fwrite(fid_bin, numel(simul_out.hs), int);
fwrite(fid_bin, simul_out.hs, double);

fprintf(fid_dat, 'M_sigma double\n');
fwrite(fid_bin, numel(simul_out.sigma), int);
fwrite(fid_bin, simul_out.sigma(1:3:end), double);
fwrite(fid_bin, simul_out.sigma(2:3:end), double);
fwrite(fid_bin, simul_out.sigma(3:3:end), double);

fprintf(fid_dat, 'M_damage double\n');
fwrite(fid_bin, numel(simul_out.damage), int);
fwrite(fid_bin, simul_out.damage, double);

fprintf(fid_dat, 'M_divergence_rate double\n');
fwrite(fid_bin, numel(simul_out.divergence_rate), int);
fwrite(fid_bin, simul_out.divergence_rate, double);

fprintf(fid_dat, 'M_h_ridged_thin_ice double\n');
fwrite(fid_bin, numel(simul_out.h_ridged_thin_ice), int);
fwrite(fid_bin, simul_out.h_ridged_thin_ice, double);

fprintf(fid_dat, 'M_h_ridged_thick_ice double\n');
fwrite(fid_bin, numel(simul_out.h_ridged_thick_ice), int);
fwrite(fid_bin, simul_out.h_ridged_thick_ice, double);

fprintf(fid_dat, 'M_random_number double\n');
fwrite(fid_bin, numel(simul_out.random_number), int);
fwrite(fid_bin, simul_out.random_number, double);

fprintf(fid_dat, 'M_Tice_0 double\n');
fwrite(fid_bin, numel(simul_out.tsurf), int);
fwrite(fid_bin, simul_out.tsurf, double);

fprintf(fid_dat, 'M_sst double\n');
fwrite(fid_bin, numel(simul_out.sst), int);
fwrite(fid_bin, simul_out.sst, double);

fprintf(fid_dat, 'M_sss double\n');
fwrite(fid_bin, numel(simul_out.sss), int);
fwrite(fid_bin, simul_out.sss, double);

fprintf(fid_dat, 'M_VT double\n');
fwrite(fid_bin, numel(simul_out.VT), int);
fwrite(fid_bin, simul_out.VT(1:2:end), double);
fwrite(fid_bin, simul_out.VT(2:2:end), double);

fprintf(fid_dat, 'M_VTM double\n');
fwrite(fid_bin, numel(simul_out.VTold), int);
fwrite(fid_bin, simul_out.VTold(1:2:end), double);
fwrite(fid_bin, simul_out.VTold(2:2:end), double);

fprintf(fid_dat, 'M_VTMM double\n');
fwrite(fid_bin, numel(simul_out.VToldold), int);
fwrite(fid_bin, simul_out.VToldold(1:2:end), double);
fwrite(fid_bin, simul_out.VToldold(2:2:end), double);

fprintf(fid_dat, 'M_UM double\n');
fwrite(fid_bin, numel(simul_out.UM), int);
fwrite(fid_bin, simul_out.UM(1:2:end), double);
fwrite(fid_bin, simul_out.UM(2:2:end), double);

fclose(fid_dat);
fclose(fid_bin);




