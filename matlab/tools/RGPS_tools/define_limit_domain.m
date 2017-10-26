function [ limit_domain ] = define_limit_domain( limit_domain_model )
%-------------- Limit the domain -------------------
% limit_domain describes a convex polygon
% or a series of convex polygons
% that will define the region where data are taken
%
% limit_domain_model    % 0  take the whole domain, 
%                       % 1 limit the domain to a polygon
%                       % 2 limit the domain to the area covered by the RGPS streams for the winter 2007-2008
%                       % 3 limit the domain to the area covered by the
%                       RGPS streams for the winter 2007-2008 and remove
%                       coastal areas


switch(limit_domain_model)
    case 0,
        limit_domain={};
    case 1,
        limit_domain(1).x=[-2050,-1700,-1000,300,700,700,500,-300,-1400,-2000,-2050];
        limit_domain(1).y=[300,1200,1500,1300,400,-200,-600,-600,-300,0,300];
    case 2,
        load('RGPS_2007-12-01_2008-06-01_traj.mat')
        
        for i=1:length(out)
            tmp_out=out{i};
            m_proj('Stereographic','lon',-45,'lat',90,'radius',60);
            
            lat=[tmp_out.meta.s_w_lat,tmp_out.meta.s_e_lat,tmp_out.meta.n_e_lat,tmp_out.meta.n_w_lat,tmp_out.meta.s_w_lat];
            lon=[tmp_out.meta.s_w_lon,tmp_out.meta.s_e_lon,tmp_out.meta.n_e_lon,tmp_out.meta.n_w_lon,tmp_out.meta.s_w_lon];
            
            %                 if i==1
            %                     figure(2)
            %                     figure(1)
            %                     worldmap([65 90], [180 360]);
            %                     geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
            %                 end
            %                 figure(1)
            %                 geoshow(lat,lon)
            
            [x,y]=m_ll2xy(lon,lat);
            x_corner_stream=(x*6378.273);
            y_corner_stream=(y*6378.273);
            
                            figure(2)
                            plot(x_corner_stream,y_corner_stream);
                            hold on
            
            limit_domain(i).x=x_corner_stream;
            limit_domain(i).y=y_corner_stream;
        end
     case 3, % a slightly smaller region to ba at least at 100 km from the coastsplots_def   
        limit_domain(1).x=[-2000,-1600,-1000,200,700,700,500,-300,-1400,-2000];
        limit_domain(1).y=[0,1100,1400,1200,400,-200,-600,-500,-200,0];
     case 4, % Rough estimate of the MY ice area
        limit_domain(1).x=[-1400-50,-1400-50,-400-50,500+50,500+50,-1400-50];
        limit_domain(1).y=[-1000-50,300+50,1100+50,1100+50,-1000-50,-1000-50];
end
%-------------- End of limit the domain -------------------

end

