function [dat] = read_SARdrift(fname)
% read the SAR drift data provided by Denis from NIERSC

m_proj('Stereographic','lon',-45,'lat',90,'radius',60);

dat.dnum = []; % end date
dat.deltat = []; % deltat in days
dat.rux = []; % velocity in m/s
dat.rvy = []; % velocity in m/s
dat.x = []; % end position
dat.y = []; %end position

showplot=0;

data=load(fname);

if(length(fname)==29)
    %   if file names include hours and minutes
    date_start=datenum(str2double(fname(1:4)),str2double(fname(5:6)),...
        str2double(fname(7:8)),str2double(fname(9:10)),...
        str2double(fname(11:12)),0);
    date_end=datenum(str2double(fname(14:17)),str2double(fname(18:19)),...
        str2double(fname(20:21)),str2double(fname(22:23)),...
        str2double(fname(24:25)),0);
elseif(length(fname)==21)
    %   If file name doesn't include hours and minutes
    date_start=datenum(str2double(fname(1:4)),str2double(fname(5:6)),...
        str2double(fname(7:8)),0,0,0);
    date_end=datenum(str2double(fname(10:13)),str2double(fname(14:15)),...
        str2double(fname(16:17)),0,0,0);
else
    disp(fname)
    error('wrong name format')
end

t0=ones(length(data),1)*date_start;
t1=ones(length(data),1)*date_end;

% LON LAT convention
[x0,y0]=m_ll2xy(data(:,1),data(:,2));
[x1,y1]=m_ll2xy(data(:,3),data(:,4));

if(min(data(:,2))<0)
    error('negative latitude: are you sure of the convention LAT LON or LON LAT')
end
    
X0=x0*6378.273;
Y0=y0*6378.273;
X1=x1*6378.273;
Y1=y1*6378.273;
Vx=(x1*6378.273-x0*6378.273)/(t1(1)-t0(1)); %velocity (x-component) in km/day
Vy=(y1*6378.273-y0*6378.273)/(t1(1)-t0(1)); %velocity (y-component) in km/day
ind=1:length(t0);

if(showplot)
    load arctic_coasts.mat;
    figure
    scatter(X1,Y1,4,Vx)
    colormap(jet(256));
    figure
    scatter(X1,Y1,4,Vy)
    colormap(jet(256));
    figure
    scatter(X1,Y1,4,(Vx.^2+Vy.^2).^(0.5))
    colormap(jet(256));
    caxis([0 40])
    colorbar
    hold on
    plot(x_coast+5,y_coast-5,'k','linewidth',0.5);
end

if (length(ind)~=length(Vx))
    warning('dimension mismatch, reducing ind to fit Vx')
    ind = 1:length(Vx)
end
    
dat.dnum = [dat.dnum;t1(ind)]; % end date
dat.deltat = [dat.deltat;t1(ind)-t0(ind)]; % deltat in days
dat.rux = [dat.rux;Vx(ind)*1000/24/60/60]; % velocity in m/s
dat.rvy = [dat.rvy;Vy(ind)*1000/24/60/60]; % velocity in m/s
dat.x = [dat.x;X1(ind)]; % end position
dat.y = [dat.y;Y1(ind)]; %end position

if(showplot)
    X1_old=X1;
    Y1_old=Y1;
    
    X1=X1(ind);
    Y1=Y1(ind);
    tri=delaunay([X1,Y1]);
    
    tr=TriRep(tri,X1,Y1);
    e=edges(tr);
    ti = edgeAttachments(tr,e);
    
    figure; clf
    length(tri)
    triplot(tr); axis image
    i_plot=0;
    hold on
    plot(X1_old,Y1_old,'or')
end

end
