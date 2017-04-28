addpath tools/subfunctions
addpath ../data
addpath ../mesh
issm  = getenv('ISSM_DIR');
addpath([issm,'/lib']);
disp('add path to starting mesh if want to plot coastlines')
