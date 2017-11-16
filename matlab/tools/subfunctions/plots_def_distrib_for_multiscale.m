function plots_def_distrib_for_multiscale(defo,invar_name,nb_bins,bin_nb)

colors={[1 0 0],[1 0.5 0],[1 1 0],[0.33 1 0],[0 1 1],[0 0 1],[0.33 0 1],[0.5 0 1],[1 0 1],[1 0.33 1]};
colors_reverse=fliplr(colors);
id_color=colors{bin_nb-1};

if(~isempty(defo.data))
    def = abs(getfield(defo.data_bin(bin_nb).invar,invar_name));
    mypdf2(def,nb_bins,{'sq-' 'color' id_color});
    %mean_def=mean(def);
    %disp([sprintf('mean %s rate: ',invar_name) num2str(mean_def)])
end
