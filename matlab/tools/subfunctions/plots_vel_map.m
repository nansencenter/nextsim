function plots_vel_map(defo_vec,title_defo,domain,masked,figure_format,date)

% We first read in the log file to know which mesh has been used
simul_in  = read_simul_in(['','nextsim.log'],0);
simul_in=simul_in.simul;
%

if(strcmp(simul_in.mesh_filename,'small_arctic_10km.msh'))
    mesh_filename='small_Arctic_10km.msh';
end
if(~exist(mesh_filename,'file'))
    mesh_filename='';
    warning(['add directory with meshfile ''',simul_in.mesh_filename,''' to path']);
end

visible = 1;
nb_defo=length(defo_vec);
if(length(title_defo)~=nb_defo)
    error('wrong title_defo')
end

for defo_nb=1:nb_defo
    defo=load_defo(defo_vec{defo_nb}); 
    
    jet_colormap=jet(256);
    
    if(isfield(defo.data,'xy_tricorner'))
        plot_tricorner(defo,'u','jet(256)',[-10 10],['u component in km/day (', title_defo{defo_nb} ,')'],domain,[],masked,visible,figure_format,date,mesh_filename);
        plot_tricorner(defo,'v','jet(256)',[-10 10],['v component in km/day (', title_defo{defo_nb} ,')'],domain,[],masked,visible,figure_format,date,mesh_filename);
        plot_tricorner(defo,'speed',jet_colormap(128:256,:),[0 20],['Speed in km/day (', title_defo{defo_nb} ,')'],domain,[],masked,visible,figure_format,date,mesh_filename);
    else
        max_speed_to_plot=10;
        
        tmp.u=defo.data.uv(:,1)/1000;
        tmp.v=defo.data.uv(:,2)/1000;
        tmp.speed=sqrt(defo.data.uv(:,1).^2+defo.data.uv(:,2).^2)/1000;
        
        if(masked==1)
            indices=defo.data.indices;
        else    
            indices=1:length(defo.data.dnum);
        end
        length(indices)
        % %
        plot_data(tmp.u(indices),defo.data.x(indices)/1000,defo.data.y(indices)/1000,title_defo{defo_nb},' ice u-velocity (km/day) from ',-max_speed_to_plot,max_speed_to_plot)
        plot_data(tmp.v(indices),defo.data.x(indices)/1000,defo.data.y(indices)/1000,title_defo{defo_nb},' ice v-velocity (km/day) from ',-max_speed_to_plot,max_speed_to_plot)
        plot_data(tmp.speed(indices),defo.data.x(indices)/1000,defo.data.y(indices)/1000,title_defo{defo_nb},' ice speed (km/day) from ',0,2*max_speed_to_plot)
    end
end

end