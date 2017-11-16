function [data_x,data_y,data_h,model_h,error,relative_error,ICESAT_campaign]=compare_ICESAT(saved_simul_in,ICESAT_campaign)

% Information on ICESAT campaigns
Laser{1}='2a';
Campaign{1}='ON03';
Period_start{1}=datenum('Sep 24 2003');
Period_end{1}=datenum('Nov 18 2003');
Days_of_Operation{1}=55;
Laser{2}='2b';
Campaign{2}='FM04';
Period_start{2}=datenum('Feb 17 2004');
Period_end{2}=datenum('Mar 21 2004');
Days_of_Operation{2}=34;
Laser{3}='3a';
Campaign{3}='ON04';
Period_start{3}=datenum('Oct 03 2004');
Period_end{3}=datenum('Nov 08 2004');
Days_of_Operation{3}=37;
Laser{4}='3b';
Campaign{4}='FM05';
Period_start{4}=datenum('Feb 17 2005');
Period_end{4}=datenum('Mar 24 2005');
Days_of_Operation{4}=36;
Laser{5}='3d';
Campaign{5}='ON05';
Period_start{5}=datenum('Oct 21 2005');
Period_end{5}=datenum('Nov 24 2005');
Days_of_Operation{5}=35;
Laser{6}='3e';
Campaign{6}='FM06';
Period_start{6}=datenum('Feb 22 2006');
Period_end{6}=datenum('Mar 27 2006');
Days_of_Operation{6}=34;
Laser{7}='3g';
Campaign{7}='ON06';
Period_start{7}=datenum('Oct 25 2006');
Period_end{7}=datenum('Nov 27 2006');
Days_of_Operation{7}=34;
Laser{8}='3h';
Campaign{8}='MA07';
Period_start{8}=datenum('Mar 12 2007');
Period_end{8}=datenum('Apr 14 2007');
Days_of_Operation{8}=34;
Laser{9}='3i';
Campaign{9}='ON07';
Period_start{9}=datenum('Oct 02 2007');
Period_end{9}=datenum('Nov 05 2007');
Days_of_Operation{9}=37;
Laser{10}='3j';
Campaign{10}='FM08';
Period_start{10}=datenum('Feb 17 2008');
Period_end{10}=datenum('Mar 21 2008');
Days_of_Operation{10}=34;

selected_campaign=0;
for i=1:length(Campaign)
    if(strcmp(lower(ICESAT_campaign),lower(Campaign{i})))
        selected_campaign=i;
        break
    end
end

if(selected_campaign)
    % load the data
    filename=['icesat_icethk_' lower(ICESAT_campaign) '_filled.dat'];
    [IceSat]=read_Icesat_data(filename);
    
    min_box_size=IceSat.x(1,2)-IceSat.x(1,1);
    xmin=IceSat.x(1,1)-min_box_size/2;
    ymin=IceSat.y(1,1)-min_box_size/2;
    xmax=IceSat.x(1,end)+min_box_size/2;
    ymax=IceSat.y(end,1)+min_box_size/2;
    ratio_dom_resol=length(IceSat.x);
    
    % load simul_in
    load(saved_simul_in);
    disp(saved_simul_in);
    
    simul_in_domain=simul_in.domain;
    simul_in_resol=simul_in.resol;
    simul_in_name=simul_in.simul_in_name;
    meshfile=[simul_in_domain simul_in_resol '.mat'];
    
    date_start=Period_start{selected_campaign};
    date_end=Period_end{selected_campaign}+1;
    
    directory_simul_out=fileparts(which(saved_simul_in));
    disp(directory_simul_out);
    
    % steps to be loaded
    output_timestep_days=simul_in.output_timestep/(24*3600);
    
    round_dnum_1=floor((date_start-simul_in.time_init)/output_timestep_days)*output_timestep_days;
    round_dnum_2=ceil((date_end  -simul_in.time_init)/output_timestep_days)*output_timestep_days;
    
    step_from = max(0,round_dnum_1/output_timestep_days);
    step_to   = round_dnum_2/output_timestep_days;
    
    % load the step 0
    saved_simul_out=['simul_out_' meshfile(1:end-4) '_' simul_in_name '_step' num2str(step_from) '.mat'];
    load(saved_simul_out);
    
    % Fill in with each step
    nb_step=step_to-step_from+1;
    
    step_ids=1:nb_step;
    
    %-------------------------------------------------------------------------%
    % Loading model data
    %-------------------------------------------------------------------------%
    model_data_nb=zeros(size(IceSat.Th));
    model_h=IceSat.Th*0;
    
    nb_steps=length(step_ids);
    for i=1:nb_steps;
        
        % load the step i
        saved_simul_out=['simul_out_' meshfile(1:end-4) '_' simul_in_name '_step' num2str(step_from+step_ids(i)-1) '.mat'];
        load(saved_simul_out);
        disp(saved_simul_out)
        
        node_x=simul_out.bamg.mesh.Vertices(:,1)*1000+simul_out.UM(1:2:end-1);
        node_y=simul_out.bamg.mesh.Vertices(:,2)*1000+simul_out.UM(2:2:end);
        
        elements_x=mean(node_x(simul_out.bamg.mesh.Triangles(:,1:3)),2)/1000;
        elements_y=mean(node_y(simul_out.bamg.mesh.Triangles(:,1:3)),2)/1000;
        
        within_domain=(elements_x>xmin).*(elements_x<xmax).*(elements_y>ymin).*(elements_y<ymax);

        cells_x=elements_x (within_domain==1); % km
        cells_y=elements_y (within_domain==1); % km
        if(isfield(simul_out,'h_thin'))
            cells_h=simul_out.h(within_domain==1)+simul_out.h_thin(within_domain==1); % m
        else
            cells_h=simul_out.h(within_domain==1); % m
        end
        indices_x=(floor((cells_x-xmin)/min_box_size)+1)';
        indices_y=(floor((cells_y-ymin)/min_box_size)+1)';
        
        indices_tot=(indices_x-1)*ratio_dom_resol+indices_y;
                
        % define a list of unique box indices
        [indice_unique,i_unique_tot,i_tot_unique] = unique(indices_tot);
        
        counter_tmp=1;
        for j=1:length(indice_unique)
            
            % if ice in the data
            data_thickness=IceSat.Th(indice_unique(j));
            if(~isnan(data_thickness) && data_thickness>=0)
                model_data_nb(indice_unique(j)) = model_data_nb(indice_unique(j))+1;
                model_h(indice_unique(j)) = model_h(indice_unique(j))+cells_h(i_unique_tot(j));
            end
        end    
    end
    
    model_h(model_data_nb>0)=model_h(model_data_nb>0)./model_data_nb(model_data_nb>0);
    
    model_h(IceSat.Th==-1)=-1/100;
    
    data_x=IceSat.x;
    data_y=IceSat.y;
    data_h=IceSat.Th/100;
    
    figure
    h=pcolor(data_x,data_y,data_h);
    caxis([0, 4]);
    set(h, 'EdgeColor', 'none');
    colorbar
    colormap(jet)
    title(['IceSAT (',ICESAT_campaign,')'])
    
    figure
    h=pcolor(data_x,data_y,model_h);
    caxis([0, 4]);
    set(h, 'EdgeColor', 'none');
    colormap(jet)
    colorbar
    title(['Simulated ice thickness (',ICESAT_campaign,')'])
    
    figure
    error=(model_h-data_h);
    h=pcolor(data_x,data_y,error);
    caxis([-2, 2]);
    set(h, 'EdgeColor', 'none');
    colormap(jet)
    colorbar
    title(['Error (',ICESAT_campaign,')']) 
    
    figure
    relative_error=(model_h-data_h)./(data_h+1e-9);
    relative_error(data_h<0)=nan;
    h=pcolor(data_x,data_y,relative_error);
    caxis([-1, 1]);
    set(h, 'EdgeColor', 'none');
    colormap(jet)
    colorbar
    title(['Relative error (',ICESAT_campaign,')']) 
else
    error(['no ICESAT campaign corresponds to ' ICESAT_campaign])
end
