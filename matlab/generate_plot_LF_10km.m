% Draw a map of observed and simulated lead fraction on 2011-03-14

mppfile = which('NpsNextsim.mpp');
m_proj('stereographic','lon',225-45,'lat',84,'rad',20, 'rot', 225)

load('ice_conc_cmap64.mat')
%% Load data
for date=datenum(2007,04,01):datenum(2007,04,14)
    close all
    i=date-datenum(2006,11,15);
    
    load(['/Volumes/SYLVAIN_EXTERNAL_MEMORY/local-data/LF_newTP_dataset/data_mat/LF_newTP_' datestr(date, 'YYYYmmdd')])
    lon_lf = lon;
    lat_lf = lat;
    
    LFHQ = LF_corr/100;
    LFHQ(quality_flag~=0) = nan;
    LFHQ(LFHQ>1) = nan;
    
    [xobs, yobs] = mapx_forward(mppfile, lon(:)', lat(:)');
    
    LFHQ_1=LFHQ(1:end-1,1:end-1);
    LFHQ_2=LFHQ(2:end,1:end-1);
    LFHQ_3=LFHQ(1:end-1,2:end);
    LFHQ_4=LFHQ(2:end,2:end);
    LFHQ_max=max([LFHQ_1(:),LFHQ_2(:),LFHQ_3(:),LFHQ_4(:)],[],2);
    LFHQ_max=reshape(LFHQ_max,size(LFHQ,1)-1,size(LFHQ,2)-1);
    LFHQ(1:end-1,1:end-1)=LFHQ_max;
    
    %subplot(1,2,1)
    % axes(ah(1))
    figure(1)
    m_pcolor(lon_lf,lat_lf,double(LFHQ>8e-3)), shading flat %, colorbar
    hold off
    colormap(ice_conc_cmap64);
    m_coast('patch',[.7 .7 .7]);
    m_grid('xticklabels',[],'yticklables',[])
    %caxis([0 0.01])
    title(['Observations ' datestr(date)])
    %
    % print('../figs/map.pdf','-dpdf','-r600')
    print(['figs_10km/LF_map_' num2str(i) '.png'],'-dpng','-r600')
end
