
function plots_def_map(defo_vec,title_defo,domain,masked,figure_format,date)

% We first read in the log file to know which mesh has been used
simul_in  = read_simul_in(['','nextsim.log'],0);
simul_in=simul_in.simul;
%
if(strcmp(simul_in.mesh_filename,'small_arctic_10km.msh'))
    mesh_filename='small_Arctic_10km.msh';
end
if ~exist(mesh_filename,'file')
    mesh_filename='';
    warning(['add directory with meshfile ''',mesh_filename,''' to path']);
end

nb_defo=length(defo_vec);
if(length(title_defo)~=nb_defo)
    error('wrong title_defo')
end

for defo_nb=1:nb_defo
    %defo=defo_vec(defo_nb);
    defo=load_defo(defo_vec{defo_nb}); 
%    plot_tricorner(defo,'eps','marsan_2004',[0 0.04],['Total deformation in /day (', title_defo{defo_nb} ,')'],domain,[],masked,1,figure_format,date,mesh_filename);

    plot_tricorner(defo,'div','blue2red',[-0.04 0.04],['Divergence rate in /day (', title_defo{defo_nb} ,')'],domain,[],masked,1,figure_format,date,mesh_filename);
    plot_tricorner(defo,'shear','gray2red',[0.005 0.1],['Shear rate in /day (', title_defo{defo_nb} ,')'],domain,[],masked,1,figure_format,date,mesh_filename);
    plot_tricorner(defo,'vor','kwok',[-0.06 0.06],['Vorticity in /day (', title_defo{defo_nb} ,')'],domain,[],masked,1,figure_format,date,mesh_filename);
end

end