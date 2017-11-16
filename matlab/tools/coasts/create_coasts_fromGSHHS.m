[ncst,Area,k]=mu_coast('l','gshhs_l.b');

x_coast=zeros(length(ncst),1);
y_coast=zeros(length(ncst),1);

angle_stereo_mesh=-45;
m_proj('Stereographic','lon',angle_stereo_mesh,'lat',90,'radius',40);

for i=1:length(ncst)
    if isnan(ncst(i,1))==0
        [a,b]=m_ll2xy(ncst(i,1),ncst(i,2));
    else
        a=NaN;b=NaN;
    end;
    x_coast(i)=a;
    y_coast(i)=b;
    100*i/length(ncst)
end;
x_coast=x_coast*6378.273;
y_coast=y_coast*6378.273;
lat_coast=ncst(:,2);
lon_coast=ncst(:,1);

clear Area a b ans i k ncst
save('arctic_coasts.mat',x_coast,y_coast,lat_coast,lon_coast)