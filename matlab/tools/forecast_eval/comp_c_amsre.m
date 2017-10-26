function [c,c_mod2amsre,x,y]=comp_c_amsre(filelist,datafile1,datafile2,varargin)
%Phil: interpolate simul_outs to AMSRE ice concentraion
%Should be give a fileslist with all the simul_out files which are all
%interpolated to the AMSRE grid and then averaged before comparing. 
%We will first load the first simulout and the satelitte data, after which
%we will do the remaining simul_outs in the filelist.
%output is the amsre concentration and the aveage simul_out concentration
%averaged to the AMSRE grid. 

%todo
%include a cutoff thickness

%current version is not finished and in beta stage

%area_box gives an x,y box to which the data is cut

%optional input
nVarargs = length(varargin);
if nVarargs >= 1, area_box = varargin{1}; end
if nVarargs >= 2, error('Too many inputs'), end

cutoff_thickness = 0.2 %m

%get the amount of files for later book keeping
[bla,nfile]=size(filelist);


disp('number of files being loaded:')
nfile


%Loading the first file just for the mesh 
load(filelist{1});
[mesh] = importbamg(simul_out.bamg.mesh, simul_out.bamg.geom);%getting the mesh





%for the AMSRE concentration:
c=double(hdfread(datafile1,'/ASI Ice Concentration', 'Index', {[1 1],[1 1],[1792 1216]}));

%charge un carte de thin ice concentration (from AMSRE, 6.5km)
try
    data=netcdf(datafile2);
    lat=double(data.VarArray(1,1).Data); %in degree north
    lon=double(data.VarArray(1,2).Data); %in degree East
    time=data.VarArray(1,3).Data; %in days since 2000.01.01
    tic=double(squeeze(data.VarArray(1,4).Data)); %in percent
    flag=double(squeeze(data.VarArray(1,5).Data)); % 0=watermask  1=ice(data ok) 100=landmask 110=ice(datanot ok) 120=missing value
    scale_factor=double(data.VarArray(1,4).AttArray(1,7).Val);
    tic=tic*scale_factor;
catch err
    if ~strcmp(err.identifier, 'MATLAB:nonStrucReference'), rethrow(err), end
    warning('No lead fraction available for this date.')
    disp('Using AMSRE__LeadFraction__UHAM-CliSAP-ICDC__v01__6.25km__20080101.nc for lon, lat, and a zero for lead fraction.')
    data=netcdf('AMSRE__LeadFraction__UHAM-CliSAP-ICDC__v01__6.25km__20080101.nc');
    lat=double(data.VarArray(1,1).Data); %in degree north
    lon=double(data.VarArray(1,2).Data); %in degree East
    % Set no lead fraction
    tic = zeros(data.DimArray(2).Dim,data.DimArray(1).Dim);
    flag = ones(data.DimArray(2).Dim,data.DimArray(1).Dim);
end


%we remove the bad quality data and apply the land-sea ice masks
f=find(flag==100  | flag==120); % | flag==110
tic(f)=NaN;
clear f;

f=find(flag==110);
tic(f)=0.;
clear f;

%we remove the bad quality data and apply the land-sea ice masks
f=find(flag==0); % | flag==110
tic(f)=100;
clear f;



%fill up the north pole
%Shouldn't be needed this way around Phil
% f=find(lat>=85);
% tmp=c(f);
% tmp2=tic(f);
% g=find(isnan(tmp)==1);
% tmp(g)=100;
% g=find(isnan(tmp2)==1);
% tmp2(g)=0;
% c(f)=tmp;
% tic(f)=tmp2;



%masking the concentration
mask=isnan(c);


%setting the projection and 
m_proj('Stereographic','lon',-45,'lat',90,'radius',60);
[x,y]=m_ll2xy(lon,lat);
x=(x*6378.273);%-20;
y=(y*6378.273);%-10;


%Removing something, not 100% sure what
f=find(x>=1900 & x<=2500);
g=find(y(f)>=-2400 & y(f)<-200);
c(f(g))=0;
tic(f(g))=0;

%using mask to revert changes done in last step? Phil
c(mask)=NaN;
tic(mask)=NaN;


%Seens to just cut off things far beyond the model grid, Phil

%lets cut the axtra 50 for shits and giggles
s=size(x);
line=ceil(s(1)/2);
%f=find(x(line,:)>=min(mesh.element.x)-50 & x(line,:)<=max(mesh.element.x)+50);
f=find(x(line,:)>=min(mesh.element.x) & x(line,:)<=max(mesh.element.x));
x=x(:,f(1):f(end));
y=y(:,f(1):f(end));
c=c(:,f(1):f(end));
tic=tic(:,f(1):f(end));

%Adding mask cut off
mask=mask(:,f(1):f(end));


s=size(y);
col=ceil(s(2)/2);
%f=find(y(:,col)>=min(mesh.element.y)-50 & y(:,col)<=max(mesh.element.y)+50);
f=find(y(:,col)>=min(mesh.element.y) & y(:,col)<=max(mesh.element.y));
y=y(f(1):f(end),:);
x=x(f(1):f(end),:);
c=c(f(1):f(end),:);
tic=tic(f(1):f(end),:);

%Adding mask cut off
mask=mask(f(1):f(end),:);



%c_mod2amsremean
c_m2a_mean = c*0.;
h_m2a_mean = c*0.;
clear simul_out mesh
   
%Now we will load all of the simul_out files, regridd them, and merge them
for i=1:nfile

   load(filelist{nfile});
   c_model=simul_out.c;
   h_model=simul_out.h;
   [mesh] = importbamg(simul_out.bamg.mesh, simul_out.bamg.geom);%getting the mesh
   
   c_mod2amsre=(griddata(mesh.element.x,mesh.element.y,c_model,x,y,'linear'));
   h_mod2amsre=(griddata(mesh.element.x,mesh.element.y,h_model,x,y,'linear'));
   
   clear simul_out mesh
   
   
   %Adding the data together to mean it
   c_m2a_mean = c_m2a_mean + c_mod2amsre;
   h_m2a_mean = h_m2a_mean + h_mod2amsre;
end


%Applying land mask of observations
c_m2a_mean(mask)=NaN;
h_m2a_mean(mask)=NaN;
   
%Calculate mean
c_m2a_mean = c_m2a_mean/real(nfile);
h_m2a_mean = h_m2a_mean/real(nfile);



%Now we cut off everything where the model thickness is between 0 and
%cutoff_thickness

thick=h_m2a_mean./c_m2a_mean;

thick_mask=find(thick>0.00001 & thick<cutoff_thickness);

plot_matrix = c*0.;
plot_matrix(thick_mask)=c_m2a_mean(thick_mask);
figure()
imshow(flipud(plot_matrix))
  
c_m2a_mean(thick_mask)=NaN;



% %interpolation from simul_out to AMSRE grid
% 
% c_mod2amsre=(griddata(mesh.element.x,mesh.element.y,c_model,x,y,'linear'));
% 
% %Applying land mask of observations
% c_mod2amsre(mask)=NaN;



%masking the NaNs of the model
mask_m2a=isnan(c_m2a_mean);

c(mask_m2a)=NaN;

%moving observations from percent to fraction
c=c/100.;





if(exist('area_box'))

    %If the are_box exists we cut both c and c_m2a_mean to it. 
    %Well, lets just turn everything outside into NANs, Phil
    
    f1=find(x< area_box(1));
    f2=find(x> area_box(2));
    f3=find(y< area_box(3));
    f4=find(y> area_box(4));
    c(f1) = NaN;
    c(f2) = NaN;
    c(f3) = NaN;
    c(f4) = NaN;
    c_m2a_mean(f1) = NaN;
    c_m2a_mean(f2) = NaN;
    c_m2a_mean(f3) = NaN;
    c_m2a_mean(f4) = NaN;
 end





%  figure()
%  imshow(flipud(c))
  
%  figure()
%  imshow(flipud(thick_mask))
  
%   figure()
%   imshow(flipud(c_m2a_mean))
%   
%   figure()
%   imshow(flipud(c_m2a_mean-c))

return 
