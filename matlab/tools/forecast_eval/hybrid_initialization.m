function [simul_out,simul_in] = hybrid_initialization(simul_in,varargin)
%*************************************************************************%
%* Merges model output with imported data to form the intial condition for
%a new forecast.
%* Called when use_simul_out ==2
%* if no past forecasts are found, use_simul_out is set to 0 and the
%traditional intitialization is done in nextsim
%*************************************************************************%

%The new state is the forecast + weight*(obs-forecast)
%weight = error.fore/(error.for+error.obs)
%Errors are absolute errors
%Errors are calculated for each cell separately.
%Errors increase with forecasts age (i.e. how far the forecast was
%integrated into the futur) by being multiplied with error.day.param
%The same would apply to observations, but for the moment we assume all
%data is 'fresh'
%For each day the observation is old the weight is reduced based on an
%estimated change rate of the variable error.day.h weight = error.day^(number of days delay)*weight

%Settings were moved to simul_in 
%use Andrea SMOS correction
%Increases 
%SMOS_correction = 1; %1 to use correction
SMOS_correction = simul_in.fore.SMOS_correction;

%increase concentration uncertainty for thin ice 
%thin_ice_uncertainty = 1; %1 to increase uncertainty
thin_ice_uncertainty = simul_in.fore.thin_ice_uncertainty; 


%Estimated error values, moved to simul_in
error.for.thick=simul_in.fore.error.for.thick; %= 0.1;
error.for.c    =simul_in.fore.error.for.c    ; %= 0.2;
error.for.d    =simul_in.fore.error.for.d    ; %= 0.2;
               
error.day.thick=simul_in.fore.error.day.thick; %= 1.1;
error.day.c    =simul_in.fore.error.day.c    ; %= 1.3;
error.day.d    =simul_in.fore.error.day.d    ; %= 0.2;
               
error.obs.thick=simul_in.fore.error.obs.thick; %= 1;
error.obs.c    =simul_in.fore.error.obs.c    ; %= 0.2;
error.obs.d    =simul_in.fore.error.obs.d    ; %= 1;



verbosity_flag =  0; %0 = basic output, 1 increased output
nVarargs = length(varargin);
if nVarargs >= 1, verbosity_flag = varargin{1}; end
if nVarargs >= 2, error('Too many inputs'), end


nb_past            = 3; %How many days in the past are looked through if the last forecast is not available
force_restart_flag = 0; %If no previous forecast can be found, fore_restart_flag will be set to 1 and the

%First step find the output file from a day before
%If there is no data from the last forecast avaiable it will look back further and further until nb_past
%For each day back the weights are reduced, as the past forecase becomes less and less reliable
previous_name = ['forecast_' datestr((simul_in.time_init-1),'YYYYmmDD')];
step_nb       = 1;
saved_simul_out=['simul_out_' simul_in.meshfile(1:end-4) '_' previous_name '_step' num2str(1) '.mat'];

for i = 1:nb_past;
    
    previous_name = ['forecast_' datestr((simul_in.time_init-i),'YYYYmmDD')];
    step_nb       = i;
    saved_simul_out=['simul_out_' simul_in.meshfile(1:end-4) '_' previous_name '_step' num2str(1) '.mat'];
    error.for.c = error.for.c*error.day.c^(1-i);
    error.for.thick = error.for.thick*error.day.thick^(1-i);
    error.for.d = error.for.d*error.day.d^(1-i);
    
    
    try
        load(saved_simul_out)
        %make sure the initial time of the simul_in agrees with the simul_out.current_time
        if simul_in.time_init == simul_out.current_time
            break
        else
            error('the loaded forecast does not have the correct date')
        end
        
    catch
        disp('looking for an earlier forecast available for initialization')
    end
    
    %If no previous forecast is available
    if i==3
        force_restart_flag            =1;
        simul_in.use_simul_out        =0;
        simul_in.spinup_duration      =1/8*86400;
        simul_out.restart_failure     =1;
        disp('forecast not available for initialization, simul_in.use_simul_out has been set to 0')
    end
end





if force_restart_flag ==0
    
    %calculate mesh
    [mesh, simul_in.ind_node_fix_bnd, simul_in.ind_node_free_bnd] = importbamg(simul_out.bamg.mesh, simul_out.bamg.geom);
    simul_out.p=0;
    fore.c      = simul_out.c;
    fore.d      = simul_out.damage;
    
    fore.thick  = zeros(size(simul_out.c));
    nonzero     = find(fore.c>0.);
    fore.thick(nonzero) = simul_out.h(nonzero)./simul_out.c(nonzero);
    
    %setting up error arrays
    
    error.for.thick = error.for.thick + zeros(mesh.Ne,1);
    error.for.c     = error.for.c     + zeros(mesh.Ne,1);
    error.for.d     = error.for.d     + zeros(mesh.Ne,1);
   
    error.obs.thick = error.obs.thick + zeros(mesh.Ne,1);
    error.obs.c     = error.obs.c     + zeros(mesh.Ne,1);
    error.obs.d     = error.obs.d     + zeros(mesh.Ne,1);
    
    %integral quantaties of last forecast for monitoring
    total_volume_pre      = simul_out.surface'*simul_out.h;
    total_snow_volume_pre = simul_out.surface'*simul_out.hs;
    total_area_pre        = simul_out.surface'*simul_out.c;
    
    %---------------------------------------------------------------------%
    % Reading in ice concentration observation                                        %
    %---------------------------------------------------------------------%
    switch simul_in.init_conc
        case 5
            % we load ice concentration from AMSR2
            fprintf('Ice concentration initialization from AMSR2 (concentration only)...\n\n');
            try
                obs.c=get_iceconcAMSR2(simul_in.time_init,mesh.element);
            catch
                fprintf('no AMSR2 data found, using previous forecast instead...\n\n');
                obs.c=simul_out.c;
            end
            
            %Satellite concntration is not reliable for thickness below 15 cm
            %We don't try to correct the bias, but we increase the concentration uncertainty for thin ice
            %Based on Natalias paper, we increase the uncertainity by sqrt(2) below 0.15 m and by 2 for ice thinner than 0.075m
            if thin_ice_uncertainty == 1; %1 to use correction
               thin_ice = find(fore.thick<0.15);
               error.obs.c(thin_ice) = error.obs.c(thin_ice).*sqrt(2);
               
               thin_ice = find(fore.thick<0.075);
               error.obs.c(thin_ice) = error.obs.c(thin_ice).*sqrt(2);
               fprintf('Ice concentration initializated and thin ice uncertainity added \n\n');
            else
               fprintf('Ice concentration initializated\n\n');    
            end

            
        case 6
            % we load ice concentration from TOPAZ4 forecast
            fprintf('Ice concentration initialization from TOPAZ4 forecast...\n\n');
            try
                obs.c=get_iceconcTP4_forecast(simul_in.time_init,mesh.element);
            catch
                fprintf('no TOPAZ data found, using previous forecast instead...\n\n');
                obs.c=simul_out.c;
            end
            
        otherwise
            error('Your choice for simul_in.init_conc and use_simul_out==2 does not correspond to a valid option')
    end
    
    %---------------------------------------------------------------------%
    % Reading in ice thickness observation                                        %
    %---------------------------------------------------------------------%
    switch simul_in.init_thick
        
        case 0 %Would normally be a constant thickness, but here we just leave the thickness untouched
            obs.thick = 1.*fore.thick;
            
        case 1
            % we load ice thickness from SMOS
            %fprintf('Ice thickness initialization from SMOS...\n\n');
            for i = 1:nb_past;
                
                try
                    if SMOS_correction == 1
                        [obs.thick, error.obs.thick] = get_icethickSMOS_v3(simul_in.time_init+1-i,mesh.element,fore.c);
                        arbitraty_smos_parameter = 1.5;
                    else
                        [obs.thick, error.obs.thick] = get_icethickSMOS_v2(simul_in.time_init+1-i,mesh.element);
                        arbitraty_smos_parameter = 2.;
                    end
                    %We do not trust the accuracy of SMOS, so we set it to be a minimum of 0.1
                    error.obs.thick(find(error.obs.thick<0.1)) = 0.1;
                    
                    
                    %We are currently dealing with the problem that the SMOS errors
                    %are for the SMOS pixel, not the smaller neXtSIM pixels.
                    %As a temporary fix we assume that the error is always at least
                    %0.1, and multiply it by a arbitrary parameter
                    error.obs.thick = error.obs.thick.*arbitraty_smos_parameter;
                    
                    %and now the added error for older SMOS values
                    error.obs.thick = error.obs.thick.*error.day.thick^(i-1);
                    
                    fprintf(['SMOS data taken from ' num2str(i-1) ' day ago ...\n\n']);
                    break
                end
                if(i==nb_past)
                    fprintf('no SMOS data found, using previous forecast instead ...\n\n');
                    obs.thick=simul_out.thick;
                end
            end
            
            
            
            
            
        case 2
            % we load ice thickness from TOPAZ4
            try
                % If simul_in.init_thick == 4 we correct volume using PIOMAS
                
                obs.thick = get_icethickTP4(simul_in.init_topaz_file,simul_in.time_init,mesh.element);
                fprintf('Ice thickness initialization from TOPAZ reanalysis...\n\n');
            catch
                fprintf('no TOPAZ data found, using previous forecast instead...\n\n');
                obs.thick=simul_out.thick;
            end
            
        case 6
            % we load ice thickness from TOPAZ4 forecast 
            try
                 obs.thick = get_icethickTP4_forecast(simul_in.time_init,mesh.element);
                 fprintf('Ice thickness initialization from TOPAZ forecast...\n\n');
            catch
                fprintf('no TOPAZ data found, using previous forecast instead...\n\n');
                obs.thick=simul_out.thick;
            end
            
            
        otherwise
            error('Your choice for simul_in.init_thick and use_simul_out =2 does not correspond to a valid option')
    end;
    
    %---------------------------------------------------------------------%
    % Initial ice damage                                                  %
    %---------------------------------------------------------------------%
    
    switch simul_in.init_damage
        case 0
            fprintf('Damage reduction...\n\n');
            obs.d=0.0*ones(size(mesh.element.x));
        case 1
            % Ice damage initialisation from RGPS
            %[simul_out]=damage_from_RGPS_deformation(simul_in,element,mesh,simul_out);
            error('Sorry, not implemented yet')
        otherwise
            error('Your choice for simul_in.init_damage and use_simul_out =2 does not correspond to a valid option')
    end
    
   
    %---------------------------------------------------------------------%
    % Initial snow depth                                                  %
    %---------------------------------------------------------------------%
    %Only needed because we currently do not have precipitation data for EC
    %forecasts
    if simul_in.use_thermo_forcing == 5 & simul_in.init_snow ==1
        % We use Warren but arbitrarily we limit the snow thickness to 20% of ice thickness
        fprintf('Snow thickness initialization from Warren climatology...\n\n');
        simul_out.hs = min(0.2*simul_out.h, Warren(simul_in.time_init,mesh.element).*simul_out.c);
    end
    
    
    
    %Computing ratio of h_ridged_thick_ice to h from the forecast
    ratio_ridged  = zeros(mesh.Ne,1);
    ratio_ridged(find(simul_out.h>0)) = simul_out.h_ridged_thick_ice(find(simul_out.h>0))./simul_out.h(find(simul_out.h>0));
    
    
    %Calculate weighting from errors
    weight.c = error.for.c./(error.for.c+error.obs.c);
    weight.d = error.for.d./(error.for.d+error.obs.d);
    weight.thick = error.for.thick./(error.for.thick+error.obs.thick);
    
    
    %Now to merge the fields according to the weighting
    thick = fore.thick + weight.thick.*(obs.thick - fore.thick);
    c     = fore.c     + weight.c.*(obs.c - fore.c);
    d     = fore.d     + weight.d.*(obs.d - fore.d);
    
    
    %Calculating h 
    simul_out.h=c.*thick ;
    
    %Assigning c and d 
    simul_out.c      = c;
    simul_out.damage = d;
    
    
    %CLEANUP
    %Set thick 0.05 where c is > zero and thick == 0 
    no_h               = find(thick==0. & c > 0.); 
    simul_out.h(no_h)  = 0.05.*c(no_h);
    %Set thick and concentration to 0.00 where c is < 0.01  
    no_c               = find(simul_out.c < 0.01); 
    simul_out.h(no_c)  = 0.;
    simul_out.c(no_c)  = 0.;
    
    
    %h_ridged is scaled with the new simul_out to ensure that h_ridged
    %doesnot exceed h
    simul_out.h_ridged_thick_ice = ratio_ridged.* simul_out.h;
    
    
    
    
    %set snow to zero where c = 0
    simul_out.hs(find(simul_out.c==0 & simul_out.hs>0.))=0;
    
    
    %integral quantaties of updated initial condition
    total_volume_post      =simul_out.surface'*simul_out.h;
    total_snow_volume_post =simul_out.surface'*simul_out.hs;
    total_area_post        =simul_out.surface'*simul_out.c;
  
    fprintf('Initialization is completed \n');
    fprintf('relative changes in ice volume, area, and snow volume during merging \n');
    disp(['volume:' num2str((total_volume_post- total_volume_pre)/total_volume_pre)])
    disp(['area:'   num2str((total_area_post  - total_area_pre)/total_area_pre)])
    disp(['snow volume:' num2str((total_snow_volume_post- total_snow_volume_pre)/total_snow_volume_pre)])
    
end

end
