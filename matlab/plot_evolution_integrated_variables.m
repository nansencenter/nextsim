function [ output_args ] = plot_evolution_integrated_variables( fields,step_start, step_end, is_sequential, dirname )
%% CALL: plot_evolution_integrated_variables( field,step_start, step_end, is_sequential, dirname )
%% example of usage:
%%    plot_evolution_integrated_variables( {'Concentration'},0, 10, true, '' )
%%
%% field:
%%    Element_area
%%    M_node_max_conc 
%%    M_VT, M_VTu, M_VTv
%%    Concentration 
%%    Thickness 
%%    Snow
%%    Damage
%%    Lambda
%%    Viscosity
%%    Freezing_Temperature
%%    Tice_0
%%    SST
%%    SSS
%%    M_wind
%%    M_tair
%%    M_mslp
%%    M_Qsw_in
%%    M_tcc 
%%    M_precip 
%%    M_dair
%%    M_ocean
%%    M_ssh
%%    M_ocean_temp
%%    M_ocean_salt
%%    M_mld
%%    M_element_depth
%%    Stresses, Stresses_x, Stresses_y
%%    Nfloes
%%    Dfloe
%%    Sigma1
%%    Sigma2
%%    mld
%%    Wind
%%    Ocean
%%    Vair_factor
%%    Voce_factor
%%    bathy
%%
%% step_start: eg if it's 0, want to start the evolutiuon from step 0
%%
%% step_end: eg if it's 10, want to end the evolutiuon at step 10
%%
%%
%% is_sequential==0: results are from MPI; 1: single processor
%%
%% dir: folder containing the outputs on the mesh
%%
%% [OPTIONAL]
%% plot_options

if ~exist('dirname','var'), dirname='.'; end

if(~isempty(dirname)&& dirname(end)~='/')
    dirname=[dirname, '/'];
end
simul_in=read_simul_in([dirname 'nextsim.log' ],0);

nb_fields=length(fields);
nb_steps=step_end-step_start+1;
integrated_fields=zeros(nb_fields,nb_steps);

for p=0:0
    for i=1:nb_steps,
      step=step_start+i-1;
      
      
      if(is_sequential)
          [mesh_out,data_out] = neXtSIM_bin_revert(dirname, [], step);
      else
          [mesh_out,data_out] = neXtSIM_bin_revert(dirname, p, step);
      end
      
      Element_area=get_and_check('Element_area',data_out,dirname,step);
      
      for j=1:nb_fields
          
          [field_tmp]=extract_field(fields{j},data_out,dirname,step);
          %integrated_fields(j,i)=integrated_fields(j,i)+sum(Element_area);
          integrated_fields(j,i)=integrated_fields(j,i)+(field_tmp'*Element_area);
          %integrated_fields(j,i)=integrated_fields(j,i)+(field_tmp'*Element_area)./sum(Element_area);
      end
      
    end
    
    % Plots
    for j=1:nb_fields
       figure
       %plot(step_start:step_end,integrated_fields(j,:));
       plot(step_start:step_end,integrated_fields(j,:)/integrated_fields(j,1));
       title(fields{j})
    end
    
end

