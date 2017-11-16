load arctic_coasts;
f=find(isnan(x_coast));
x=[];y=[];lat=[];lon=[];
for i=1:length(f)-1
    i
    a=f(i);
    b=f(i+1);
    if (b-a)>10
        x=[x;x_coast(a:b)];
        y=[y;y_coast(a:b)];
        lat=[lat;lat_coast(a:b)];
        lon=[lon;lon_coast(a:b)];
    end;
end;
x_coast=x;
y_coast=y;
lat_coast=lat;
lon_coast=lon;

clear lat lon x y f i a b
save('arctic_coasts_light_tmp.mat')
    
% arctic_coast_light_tmp.mat has been modified by hand to remove rivers and some
% straits to become arctic_coast_light.mat (Ask Pierre Rampal for more information)