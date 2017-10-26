function [mesh,element]=grid2mat(filename,showplot,ocean_grid,refinement,reduced_domain,split_factor)
% Read regional.grid.a and regional.depth.a and create mesh and element
% by splitting the element of the grid into split_factor triangles (2, 4 or 8)

% Parameters that can be changed by the users

if(strcmp(ocean_grid,'TOPAZ'))
    ocean=1;
    
    if refinement
        % definition of the box within which the dist_to_coast_8 and dist_to_coast_4 is used
        % in a transition around the boxes only dist_to_coast_4 is checked
        % boxes may be overlaping
        if(split_factor~=2)
            error('refinement only works with split_factor=2')
        end
        
        % box(i,:)=[min_lat,max_lat,min_lon,max_lon];
        box(1,:)=[75,85,-130,-30];  % CAA
        box(2,:)=[67,85,-130,-70];  % CAA
        %box(3,:)=[65,75,50,60];     % kara gate
        
        % The width of the transition zone is defined in degree
        transition_width=3; % in lat or lon
        min_dist_coast_4=8; % in terms of grid cell (if there is a land cell within [i,j +- min_dist])
        min_dist_coast_8=6; %
    else
        box(1,:)=[0,0,0,0];
        transition_width=0;
        min_dist_coast_4=0; % in terms of grid cell (if there is a land cell within [i,j +- min_dist])
        min_dist_coast_8=0; % For the moment split_8 does not work because the
    end      
    
elseif(strcmp(ocean_grid,'MITgcm_9km'))
    ocean=2;
    
    if refinement
        % definition of the box within which the dist_to_coast_8 and dist_to_coast_4 is used
        % in a transition around the boxes only dist_to_coast_4 is checked
        % boxes may be overlaping
        if(split_factor~=2)
            error('refinement only works with split_factor=2')
        end
        
        % box(i,:)=[min_lat,max_lat,min_lon,max_lon];
        box(1,:)=[65,75,50,60];     % kara gate
        box(2,:)=[75,85,-130,-30];  % CAA
        box(3,:)=[60,85,-130,-80];  % CAA

        % The width of the transition zone is defined in degree
        transition_width=3; % in lat or lon
        min_dist_coast_4=12; % in terms of grid cell (if there is a land cell within [i,j +- min_dist])
        min_dist_coast_8=6; %
    else  
        box(1,:)=[0,0,0,0];
        transition_width=0;
        min_dist_coast_4=0; % in terms of grid cell (if there is a land cell within [i,j +- min_dist])
        min_dist_coast_8=0; % For the moment split_8 does not work because the
    end
  
elseif(strcmp(ocean_grid,'MITgcm_4km'))
    ocean=3;
    
    if refinement
        % definition of the box within which the dist_to_coast_8 and dist_to_coast_4 is used
        % in a transition around the boxes only dist_to_coast_4 is checked
        % boxes may be overlaping
        if(split_factor~=2)
            error('refinement only works with split_factor=2')
        end
        
        % box(i,:)=[min_lat,max_lat,min_lon,max_lon];
        box(1,:)=[65,75,50,60];     % kara gate
        box(2,:)=[75,85,-130,-30];  % CAA
        box(3,:)=[60,85,-130,-80];  % CAA

        % The width of the transition zone is defined in degree
        transition_width=3; % in lat or lon
        min_dist_coast_4=12; % in terms of grid cell (if there is a land cell within [i,j +- min_dist])
        min_dist_coast_8=6; %
    else  
        box(1,:)=[0,0,0,0];
        transition_width=0;
        min_dist_coast_4=0; % in terms of grid cell (if there is a land cell within [i,j +- min_dist])
        min_dist_coast_8=0; % For the moment split_8 does not work because the
    end
    
elseif(strcmp(ocean_grid,'simplesquare'))
    ocean=4;
    
    refinement=0;
    box(1,:)=[0,0,0,0];
    transition_width=0;
    min_dist_coast_4=0; % in terms of grid cell (if there is a land cell within [i,j +- min_dist])
    min_dist_coast_8=0; % For the moment split_8 does not work because the   
    

elseif(strcmp(ocean_grid,'archbox'))
    ocean=5;
    
    refinement=0;
    box(1,:)=[0,0,0,0];
    transition_width=0;
    min_dist_coast_4=0; % in terms of grid cell (if there is a land cell within [i,j +- min_dist])
    min_dist_coast_8=0; % For the moment split_8 does not work because the   


else
    error(['grid2mat does not support this ocean_grid:', ocean_grid])
end

if(min_dist_coast_4<=min_dist_coast_8)
   disp('min_dist_coast_4<=min_dist_coast_8, there will be no transition from split_2 to split_8') 
end

% Loading the files
switch ocean
    case 1
        % open the files
        try
            fid=fopen('regional.grid.b');
            A=fgetl(fid); 
            idm=sscanf(A,'%f');
            A=fgetl(fid); 
            jdm=sscanf(A,'%f');
            fclose(fid)
        catch
            files=dir('mask*.uf');
            error('can not get grid size');
            return
        end
        
        % we load the lat lon coordinates
        % of the center              (p point)
        % of the lower-left corner   (q point)
        % of the center of the left  (u point)
        % of the center of the lower (v point)
        try
            plon=loada('regional.grid.a',1,idm,jdm);
            plat=loada('regional.grid.a',2,idm,jdm);
            qlon=loada('regional.grid.a',3,idm,jdm);
            qlat=loada('regional.grid.a',4,idm,jdm);
            ulon=loada('regional.grid.a',5,idm,jdm);
            ulat=loada('regional.grid.a',6,idm,jdm);
            vlon=loada('regional.grid.a',7,idm,jdm);
            vlat=loada('regional.grid.a',8,idm,jdm);
        catch
            error('regional.grid.a not found');
        end
        
        try
            [pdepth]=loada('regional.depth.a',1,idm,jdm);
            pmask=double(pdepth~=pdepth(1,1));
        catch
            error('regional.depth.a not found');
        end
        
        if(reduced_domain)
            %We select the domain that will be meshed
            min_i=200;
            max_i=min(630,size(plat,1));
            min_j=435;
            max_j=size(plat,2);
        else
            % We select the whole domain
            min_i=1;
            max_i=size(plat,1);
            min_j=1;
            max_j=size(plat,2);
        end
        
        % Buldozer operations
        %We remove the few remaining grid cells lying the Mediterranean Sea
        f=find(qlat<=50 & qlon>=0 );
        pmask(f)=0;
        
        % Bering Strait is closed in TOPAZ4. To have the opportunity to open it,
        % the mask of the cells beyond Bering Strait are set to 2
        pmask(222,813:836)=2;
    case 2
        try
            plon=readbin('9km_XC.data',[420*2 384*2]);
            plat=readbin('9km_YC.data',[420*2 384*2]);
            qlon=readbin('9km_XG.data',[420*2 384*2]);
            qlat=readbin('9km_YG.data',[420*2 384*2]);
            
            % for ulon, ulat, vlon, vlat :
            m_proj('Stereographic','lon',-45,'lat',90,'radius',60);
            [qx,qy]=m_ll2xy(qlon,qlat);
            
            ux=qx;
            uy=qy;
            vx=qx;
            vy=qy;
            
            ux(:,1:end-1)=(qx(:,1:end-1)+qx(:,2:end))/2;
            uy(:,1:end-1)=(qy(:,1:end-1)+qy(:,2:end))/2;
            vx(1:end-1,:)=(qx(1:end-1,:)+qx(2:end,:))/2;
            vy(1:end-1,:)=(qy(1:end-1,:)+qy(2:end,:))/2;
            
            [ulon,ulat]=m_xy2ll(ux,uy);
            [vlon,vlat]=m_xy2ll(vx,vy);
        catch
            disp('9km_XG.data and 9km_YG.data not found -- grid files missing');
            return
        end
        
        try
            pmask =readbin('9km_hFacC.data',[420*2 384*2]);
            pdepth=readbin('9km_Depth.data',[420*2 384*2]);
        catch
            disp('9km_hFacC.data or 9km_Depth.data not found -- The mask file is missing - Operation aborted');
            return
        end
        
        if(reduced_domain)
            % We select the domain that will be meshed
            min_i=240;
            max_i=min(470,size(plat,1));
            min_j=1;
            max_j=min(200,size(plat,2));
        else
            % We select the whole domain
            min_i=1;
            max_i=size(plat,1);
            min_j=1;
            max_j=size(plat,2);
        end
    case 3
        try
            plon=readbin('XC.data',[420*4 384*4]);
            plat=readbin('YC.data',[420*4 384*4]);
            qlon=readbin('XG.data',[420*4 384*4]);
            qlat=readbin('YG.data',[420*4 384*4]);
            
            % for ulon, ulat, vlon, vlat :
            m_proj('Stereographic','lon',-45,'lat',90,'radius',60);
            [qx,qy]=m_ll2xy(qlon,qlat);
            
            ux=qx;
            uy=qy;
            vx=qx;
            vy=qy;
            
            ux(:,1:end-1)=(qx(:,1:end-1)+qx(:,2:end))/2;
            uy(:,1:end-1)=(qy(:,1:end-1)+qy(:,2:end))/2;
            vx(1:end-1,:)=(qx(1:end-1,:)+qx(2:end,:))/2;
            vy(1:end-1,:)=(qy(1:end-1,:)+qy(2:end,:))/2;
            
            [ulon,ulat]=m_xy2ll(ux,uy);
            [vlon,vlat]=m_xy2ll(vx,vy);
        catch
            disp('XG.data and YG.data not found -- grid files missing');
            return
        end
        
        try
            pmask =readbin('hFacC.data',[420*4 384*4]);
            pdepth=readbin('Depth.data',[420*4 384*4]);
        catch
            disp('hFacC.data or Depth.data not found -- The mask file is missing - Operation aborted');
            return
        end
        
        if(reduced_domain)
            % We select the domain that will be meshed
            % Phil: Reduced domain is the neXtSIMF forecast domain
            min_i=480;
            min_i=320;
            max_i=min(940,size(plat,1));
            min_j=1;
            max_j=min(400,size(plat,2));
        else
            % We select the whole domain
            min_i=1;
            max_i=size(plat,1);
            min_j=1;
            max_j=size(plat,2);
        end
        
    case 4
        m_proj('Stereographic','lon',-45,'lat',90,'radius',60);

        dx   =10 /6378.273;
        x_min=100/6378.273;
        y_min=100/6378.273;
        x_max=500/6378.273;
        y_max=500/6378.273;
                
        [px,py]=meshgrid(x_min+dx/2:dx:x_max-dx/2, y_min+dx/2:dx:y_max-dx/2)
        [qx,qy]=meshgrid(x_min     :dx:x_max-dx  , y_min     :dx:y_max-dx  )

        ux=qx;
        uy=qy;
        vx=qx;
        vy=qy;
        
        ux(:,1:end-1)=(qx(:,1:end-1)+qx(:,2:end))/2;
        uy(:,1:end-1)=(qy(:,1:end-1)+qy(:,2:end))/2;
        vx(1:end-1,:)=(qx(1:end-1,:)+qx(2:end,:))/2;
        vy(1:end-1,:)=(qy(1:end-1,:)+qy(2:end,:))/2;
        
        [plon,plat]=m_xy2ll(px,py)
        [qlon,qlat]=m_xy2ll(qx,qy)
        [ulon,ulat]=m_xy2ll(ux,uy);
        [vlon,vlat]=m_xy2ll(vx,vy);
        
        pmask=ones(size(px));
        pdepth=pmask*200; % The depth is fixed to 200 m but it could be changed 
        
        %Closed boundaries at the bottom
        pmask(1,:)                         =0;
        
        if(reduced_domain)
            % We select the domain that will be meshed
            min_i=1;
            max_i=min(10,size(plat,1));
            min_j=1;
            max_j=min(10,size(plat,2));
        else
            % We select the whole domain
            min_i=1;
            max_i=size(plat,1);
            min_j=1;
            max_j=size(plat,2);
        end
        
    case 5 %Setting up arching domain
        
        m_proj('Stereographic','lon',-45,'lat',90,'radius',60);

        dx   =5/6378.273; %Reduced later on
        x_min=100/6378.273 -dx;
        y_min=100/6378.273 -dx;
        x_max=250/6378.273 +dx;
        y_max=300/6378.273 +dx;
        
        width = 75/6378.273 ;% Dictates width of hole
        hight = 50/6378.273 ;% Dictates depth of hole
        thick = 10/6378.273 ;% Dictates depth of hole
        
        width = round(width/(2*dx))*2*dx; %Needs to be a full multiplier of dx
        hight = round(hight/dx)*dx; %Needs to be a full multiplier of dx
        thick = round(thick/dx)*dx; %Needs to be a full multiplier of dx
        [px,py]=meshgrid(x_min+dx/2:dx:x_max-dx/2, y_min+dx/2:dx:y_max-dx/2);
        [qx,qy]=meshgrid(x_min     :dx:x_max-dx  , y_min     :dx:y_max-dx  );

        ux=qx;
        uy=qy;
        vx=qx;
        vy=qy;
        
        ux(:,1:end-1)=(qx(:,1:end-1)+qx(:,2:end))/2;
        uy(:,1:end-1)=(qy(:,1:end-1)+qy(:,2:end))/2;
        vx(1:end-1,:)=(qx(1:end-1,:)+qx(2:end,:))/2;
        vy(1:end-1,:)=(qy(1:end-1,:)+qy(2:end,:))/2;
        
        [plon,plat]=m_xy2ll(px,py);
        [qlon,qlat]=m_xy2ll(qx,qy);
        [ulon,ulat]=m_xy2ll(ux,uy);
        [vlon,vlat]=m_xy2ll(vx,vy);
        
        pmask=ones(size(px));
        num_width  = round(width/dx);
        num_x      = round((x_max-x_min)/dx);
        num_y      = round((y_max-y_min)/dx);
        hight_y    = round(hight/dx)+1;
        hole_begin = num_x/2 - num_width/2 +1; 
        hole_end   = num_x/2 + num_width/2 ; 
        thick_begin= hight_y;
        thick_end  = hight_y-round(thick/dx);
        
        %Closed boundaries
        pmask(1,:)                         =0;
        pmask(:,1)                         =0;
        pmask(:,num_x)                     =0;
        pmask(num_y,:)                     =0;
        
        %Open outflow
        pmask(1,hole_begin:hole_end)       =1;
        
        %Make "land"
        %Closed case
        %pmask(1:hight_y,1:hole_begin-1)    =0;
        %pmask(1:hight_y,hole_end+1:num_x)  =0;
        
        
        %Open bottom
        pmask(hight_y,1:hole_begin-1)    =0;
        pmask(hight_y,hole_end+1:num_x)  =0;
        
        %Open whole lower domain
        pmask(1,:)       =1;
        
        pdepth=pmask*200; % The depth is fixed to 200 m but it could be changed 
        
        if(reduced_domain)
            % We select the domain that will be meshed
            min_i=1;
            max_i=min(10,size(plat,1));
            min_j=1;
            max_j=min(10,size(plat,2));
        else
            % We select the whole domain
            min_i=1;
            max_i=size(plat,1);
            min_j=1;
            max_j=size(plat,2);
        end
end

platr=plat(min_i:max_i,min_j:max_j);
plonr=plon(min_i:max_i,min_j:max_j);
qlatr=qlat(min_i:max_i,min_j:max_j);
qlonr=qlon(min_i:max_i,min_j:max_j);
ulatr=ulat(min_i:max_i,min_j:max_j);
ulonr=ulon(min_i:max_i,min_j:max_j);
vlatr=vlat(min_i:max_i,min_j:max_j);
vlonr=vlon(min_i:max_i,min_j:max_j);
maskr=pmask(min_i:max_i,min_j:max_j);
depthr=pdepth(min_i:max_i,min_j:max_j);

% size of the selected domain
[nb_x,nb_y]=size(maskr);

mask_init=((maskr~=0).*(maskr~=2));

mask_dist_4=ones(nb_x+2*min_dist_coast_4,nb_y+2*min_dist_coast_4);
mask_dist_8=ones(nb_x+2*min_dist_coast_8,nb_y+2*min_dist_coast_8);

for i=-min_dist_coast_4:min_dist_coast_4,
    for j=-min_dist_coast_4:min_dist_coast_4,
        mask_dist_4_tmp=mask_dist_4(min_dist_coast_4+1+i:min_dist_coast_4+nb_x+i,min_dist_coast_4+1+j:min_dist_coast_4+nb_y+j);
        mask_dist_4(min_dist_coast_4+1+i:min_dist_coast_4+nb_x+i,min_dist_coast_4+1+j:min_dist_coast_4+nb_y+j)=mask_dist_4_tmp.*mask_init;
    end
end

for i=-min_dist_coast_8:min_dist_coast_8,
    for j=-min_dist_coast_8:min_dist_coast_8,
        mask_dist_8_tmp=mask_dist_8(min_dist_coast_8+1+i:min_dist_coast_8+nb_x+i,min_dist_coast_8+1+j:min_dist_coast_8+nb_y+j);
        mask_dist_8(min_dist_coast_8+1+i:min_dist_coast_8+nb_x+i,min_dist_coast_8+1+j:min_dist_coast_8+nb_y+j)=mask_dist_8_tmp.*mask_init;
    end
end

% The cells at the borders of the selected domain are not meshed
% but are used to define the boundary conditions on the external borders
maskr_inner=maskr(2:end-1,2:end-1);
depthr_inner=depthr(2:end-1,2:end-1);
platr_inner=platr(2:end-1,2:end-1);
plonr_inner=plonr(2:end-1,2:end-1);
[nb_x_inner,nb_y_inner]=size(maskr_inner);
    
mask_dist_4_inner=mask_dist_4(min_dist_coast_4+2:end-min_dist_coast_4-1,min_dist_coast_4+2:end-min_dist_coast_4-1);
mask_dist_8_inner=mask_dist_8(min_dist_coast_8+2:end-min_dist_coast_8-1,min_dist_coast_8+2:end-min_dist_coast_8-1);

in_box=zeros(size(maskr_inner));
in_box_extended=zeros(size(maskr_inner));
for i=1:size(box,1)
    min_lat=box(i,1);
    max_lat=box(i,2);
    min_lon=box(i,3);
    max_lon=box(i,4);
    min_lat_extended=min_lat-transition_width;
    max_lat_extended=max_lat+transition_width;
    min_lon_extended=min_lon-transition_width;
    max_lon_extended=max_lon+transition_width;
    in_box=in_box+((platr_inner>=min_lat).*(platr_inner<=max_lat).*(plonr_inner>=min_lon).*(plonr_inner<=max_lon));
    in_box_extended=in_box_extended+((platr_inner>=min_lat_extended).*(platr_inner<=max_lat_extended).*(plonr_inner>=min_lon_extended).*(plonr_inner<=max_lon_extended));
end

if(refinement)
    % Indices of the cell to split in the inner domain
    split_cell=2*((maskr_inner~=0).*(maskr_inner~=2).*((in_box_extended==0)+((mask_dist_4_inner==1).*(mask_dist_8_inner==1)))>0);
    split_cell=split_cell+4*((maskr_inner~=0).*(maskr_inner~=2).*(mask_dist_4_inner==0).*(in_box_extended~=0).*((in_box==0)+(mask_dist_8_inner==1))>0);
    split_cell=split_cell+8*((maskr_inner~=0).*(maskr_inner~=2).*(mask_dist_8_inner==0).*(in_box~=0)>0);
    [i_2,j_2]=find(split_cell==2);
    [i_4,j_4]=find(split_cell==4);
    [i_8,j_8]=find(split_cell==8);
else
    i_2=[];
    j_2=[];
    i_4=[];
    j_4=[];
    i_8=[];
    j_8=[];
    switch split_factor
        case 2
            [i_2,j_2]=find((maskr_inner~=0).*(maskr_inner~=2));
        case 4
            [i_4,j_4]=find((maskr_inner~=0).*(maskr_inner~=2));
        case 8
            [i_8,j_8]=find((maskr_inner~=0).*(maskr_inner~=2));
        otherwise
            error(['split_factor= ' num2str(split_factor) ' is not managed by grid2mat'])
    end
end
    
if(showplot)
    figure
    [i,j]=find((pmask~=0).*(pmask~=2));
    plot(i,j,'.'); 
    hold on
    [i,j]=find((pmask==2));
    plot(i,j,'xr')
    
    figure
    plot(i_2,j_2,'.b','MarkerSize',4); 
    hold on
    plot(i_4,j_4,'.g','MarkerSize',4)
    plot(i_8,j_8,'.r','MarkerSize',4)
end

% Corresponding indices for the reduced domain 
i_2r=i_2+1;
j_2r=j_2+1;
i_4r=i_4+1;
j_4r=j_4+1;
i_8r=i_8+1;
j_8r=j_8+1;

Ncell_2=length(i_2);
Ncell_4=length(i_4);
Ncell_8=length(i_8);

% number of elements (triangles)
end_2=        Ncell_2*2;
end_4=end_2+Ncell_4*4;
end_8=end_4+Ncell_8*8;

start_2=1;
start_4=end_2+1;
start_8=end_4+1;

element.num_node=zeros(end_8,3);
element.depth=zeros(end_8,1);

if(Ncell_2>0)   
        corner_1=(i_2r  )+nb_x*(j_2r  -1);
        corner_2=(i_2r+1)+nb_x*(j_2r  -1);
        corner_3=(i_2r+1)+nb_x*(j_2r+1-1);
        corner_4=(i_2r  )+nb_x*(j_2r+1-1);
        
        random_split=rand(Ncell_2,1)>0.5;
        
        element.num_node(start_2:2:end_2-1,1)=corner_1;
        element.num_node(start_2:2:end_2-1,2)=corner_2;
        element.num_node(start_2:2:end_2-1,3)=random_split.*corner_3+(1-random_split).*corner_4;
        element.depth   (start_2:2:end_2-1,1)=depthr(corner_1);

        element.num_node(start_2+1:2:end_2,1)=corner_3;
        element.num_node(start_2+1:2:end_2,2)=corner_4;
        element.num_node(start_2+1:2:end_2,3)=random_split.*corner_1+(1-random_split).*corner_2;
        element.depth   (start_2+1:2:end_2,1)=depthr(corner_1);
end
if(Ncell_4>0)
        corner_1=          (i_4r  )+nb_x*(j_4r  -1);
        corner_2=          (i_4r+1)+nb_x*(j_4r  -1);
        corner_3=          (i_4r+1)+nb_x*(j_4r+1-1);
        corner_4=          (i_4r  )+nb_x*(j_4r+1-1);
        center  =nb_x*nb_y+(i_4r  )+nb_x*(j_4r  -1);
        
        element.num_node(start_4  :4:end_4-3,1)=corner_1;
        element.num_node(start_4  :4:end_4-3,2)=corner_2;
        element.num_node(start_4  :4:end_4-3,3)=center;
        element.depth   (start_4  :4:end_4-3,1)=depthr(corner_1);
        
        element.num_node(start_4+1:4:end_4-2,1)=corner_2;
        element.num_node(start_4+1:4:end_4-2,2)=corner_3;
        element.num_node(start_4+1:4:end_4-2,3)=center;
        element.depth   (start_4+1:4:end_4-2,1)=depthr(corner_1);
       
        element.num_node(start_4+2:4:end_4-1,1)=corner_3;
        element.num_node(start_4+2:4:end_4-1,2)=corner_4;
        element.num_node(start_4+2:4:end_4-1,3)=center;
        element.depth   (start_4+2:4:end_4-1,1)=depthr(corner_1);
        
        element.num_node(start_4+3:4:end_4  ,1)=corner_4;
        element.num_node(start_4+3:4:end_4  ,2)=corner_1;
        element.num_node(start_4+3:4:end_4  ,3)=center;
        element.depth   (start_4+3:4:end_4  ,1)=depthr(corner_1);

end
if(Ncell_8>0)
    corner_1 =            (i_8r  )+nb_x*(j_8r  -1);
    corner_2 =            (i_8r+1)+nb_x*(j_8r  -1);
    corner_3 =            (i_8r+1)+nb_x*(j_8r+1-1);
    corner_4 =            (i_8r  )+nb_x*(j_8r+1-1);
    center   =  nb_x*nb_y+(i_8r  )+nb_x*(j_8r  -1);
    mid_down =2*nb_x*nb_y+(i_8r  )+nb_x*(j_8r  -1);
    mid_up   =2*nb_x*nb_y+(i_8r  )+nb_x*(j_8r+1-1);
    mid_left =3*nb_x*nb_y+(i_8r  )+nb_x*(j_8r  -1);
    mid_right=3*nb_x*nb_y+(i_8r+1)+nb_x*(j_8r  -1);
    
    for sub_square=1:4,
        
        random_split=rand(Ncell_8,1)>0.5;
        switch sub_square
            case 1
                sub_corner=[corner_1,mid_down,center,mid_left];
            case 2
                sub_corner=[mid_down,corner_2,mid_right,center];
            case 3
                sub_corner=[center,mid_right,corner_3,mid_up];
            case 4
                sub_corner=[mid_left,center,mid_up,corner_4];
        end
        
        local_start=(sub_square-1)*2+1;
        tmp_start=start_8-1+local_start;
        tmp_end  =end_8  -8+local_start;
        element.num_node(tmp_start:8:tmp_end,1)=sub_corner(:,1);
        element.num_node(tmp_start:8:tmp_end,2)=sub_corner(:,2);
        element.num_node(tmp_start:8:tmp_end,3)=random_split.*sub_corner(:,3)+(1-random_split).*sub_corner(:,4);
        element.depth(tmp_start:8:tmp_end,1)=depthr(corner_1);
        
        local_start=(sub_square-1)*2+2;
        tmp_start=start_8-1+local_start;
        tmp_end  =end_8  -8+local_start;
        element.num_node(tmp_start:8:tmp_end,1)=sub_corner(:,3);
        element.num_node(tmp_start:8:tmp_end,2)=sub_corner(:,4);
        element.num_node(tmp_start:8:tmp_end,3)=random_split.*sub_corner(:,1)+(1-random_split).*sub_corner(:,2);
        element.depth(tmp_start:8:tmp_end,1)=depthr(corner_1);
    end
end

%---------- Build the boundary conditions

flag_fix  = 10000;
flag_free = 10001;
flag_optional = 10002;

% detection of the quadrants where the element should not be splited in 8 
% we determine the condition of the botton, right, upper and left edge
% 0 if no boundary, 1 if fix boundary, 2 if free boundary

local_start=0;

max_nb_boundary=Ncell_2*4+Ncell_4*4+Ncell_8*8;
mesh.boundary.from_msh=zeros(max_nb_boundary,3);

for i_split=[2,4,8],
    
    switch i_split
        case 2
            Ncell=Ncell_2;
            i_tmp=i_2r;
            j_tmp=j_2r;
        case 4
            Ncell=Ncell_4;
            i_tmp=i_4r;
            j_tmp=j_4r;
        case 8
            Ncell=Ncell_8;
            i_tmp=i_8r;
            j_tmp=j_8r;
    end
    
    if(Ncell>0)
        boundary_tmp=zeros(Ncell,4);
        
        i_tmp_m1=i_tmp-1;
        j_tmp_m1=j_tmp-1;
        i_tmp_p1=i_tmp+1;
        j_tmp_p1=j_tmp+1;
        
        boundary_tmp(:,1)=flag_fix*(maskr((i_tmp)+nb_x*(j_tmp_m1-1))==0);
        boundary_tmp(:,2)=flag_fix*(maskr((i_tmp_p1)+nb_x*(j_tmp-1))==0);
        boundary_tmp(:,3)=flag_fix*(maskr((i_tmp)+nb_x*(j_tmp_p1-1))==0);
        boundary_tmp(:,4)=flag_fix*(maskr((i_tmp_m1)+nb_x*(j_tmp-1))==0);
        
        boundary_tmp(:,1)=boundary_tmp(:,1)+flag_free*(maskr((i_tmp)+nb_x*(j_tmp_m1-1))~=0).*((j_tmp_m1)==1   );
        boundary_tmp(:,2)=boundary_tmp(:,2)+flag_free*(maskr((i_tmp_p1)+nb_x*(j_tmp-1))~=0).*((i_tmp_p1)==nb_x);
        boundary_tmp(:,3)=boundary_tmp(:,3)+flag_free*(maskr((i_tmp)+nb_x*(j_tmp_p1-1))~=0).*((j_tmp_p1)==nb_y);
        boundary_tmp(:,4)=boundary_tmp(:,4)+flag_free*(maskr((i_tmp_m1)+nb_x*(j_tmp-1))~=0).*((i_tmp_m1)==1   );
        
        boundary_tmp(:,1)=boundary_tmp(:,1)+flag_optional*(maskr((i_tmp)+nb_x*(j_tmp_m1-1))==2);
        boundary_tmp(:,2)=boundary_tmp(:,2)+flag_optional*(maskr((i_tmp_p1)+nb_x*(j_tmp-1))==2);
        boundary_tmp(:,3)=boundary_tmp(:,3)+flag_optional*(maskr((i_tmp)+nb_x*(j_tmp_p1-1))==2);
        boundary_tmp(:,4)=boundary_tmp(:,4)+flag_optional*(maskr((i_tmp_m1)+nb_x*(j_tmp-1))==2);
        
        % Loop to redefine the triangles in the selected quadrants
        Cell_to_visit=find(sum(boundary_tmp,2)>0);
        
        corner_1 =            (i_tmp(Cell_to_visit)  )+nb_x*(j_tmp(Cell_to_visit)  -1);
        corner_2 =            (i_tmp(Cell_to_visit)+1)+nb_x*(j_tmp(Cell_to_visit)  -1);
        corner_3 =            (i_tmp(Cell_to_visit)+1)+nb_x*(j_tmp(Cell_to_visit)+1-1);
        corner_4 =            (i_tmp(Cell_to_visit)  )+nb_x*(j_tmp(Cell_to_visit)+1-1);
        center   =  nb_x*nb_y+(i_tmp(Cell_to_visit)  )+nb_x*(j_tmp(Cell_to_visit)  -1);
        mid_down =2*nb_x*nb_y+(i_tmp(Cell_to_visit)  )+nb_x*(j_tmp(Cell_to_visit)  -1);
        mid_up   =2*nb_x*nb_y+(i_tmp(Cell_to_visit)  )+nb_x*(j_tmp(Cell_to_visit)+1-1);
        mid_left =3*nb_x*nb_y+(i_tmp(Cell_to_visit)  )+nb_x*(j_tmp(Cell_to_visit)  -1);
        mid_right=3*nb_x*nb_y+(i_tmp(Cell_to_visit)+1)+nb_x*(j_tmp(Cell_to_visit)  -1);
        
        for i=1:length(Cell_to_visit)
            
            for quadrant=1:4,
                switch quadrant
                    case 1
                        sub_corner=[corner_1(i),mid_down(i),corner_2(i)];
                    case 2
                        sub_corner=[corner_2(i),mid_right(i),corner_3(i)];
                    case 3
                        sub_corner=[corner_3(i),mid_up(i),corner_4(i)];
                    case 4
                        sub_corner=[corner_4(i),mid_left(i),corner_1(i)];
                end
                
                if(i_split==2 || i_split==4)
                    jump_corner=2;
                else
                    jump_corner=1;
                end
                
                flag_tmp=boundary_tmp(Cell_to_visit(i),quadrant);
                if(flag_tmp)
                    for j=1:jump_corner:length(sub_corner)-1;
                        local_start=local_start+1;
                        mesh.boundary.from_msh(local_start,1)=sub_corner(j);
                        mesh.boundary.from_msh(local_start,2)=sub_corner(j+jump_corner);
                        mesh.boundary.from_msh(local_start,3)=flag_tmp;
                    end
                end
            end
        end
    end
end
mesh.boundary.from_msh=mesh.boundary.from_msh(1:local_start,:);


%---------- Correct the transition between split_8 and split_4 or split_2
% Some elements from split_8 have to be merged to have a consistent transition
% from split_8 to split_2 or split_4

if(refinement)
    
    % detection of the quadrants where the element should not be splited in 8
    Quadrant_8_no_split=zeros(Ncell_8,4);
    
    i_8_m1=max(1,i_8-1);
    j_8_m1=max(1,j_8-1);
    i_8_p1=min(nb_x_inner,i_8+1);
    j_8_p1=min(nb_y_inner,j_8+1);
    
    for i=1:Ncell_8
        Quadrant_8_no_split(i,1)=(split_cell(i_8(i),j_8_m1(i))==2)||(split_cell(i_8(i),j_8_m1(i))==4);
        Quadrant_8_no_split(i,2)=(split_cell(i_8_p1(i),j_8(i))==2)||(split_cell(i_8_p1(i),j_8(i))==4);
        Quadrant_8_no_split(i,3)=(split_cell(i_8(i),j_8_p1(i))==2)||(split_cell(i_8(i),j_8_p1(i))==4);
        Quadrant_8_no_split(i,4)=(split_cell(i_8_m1(i),j_8(i))==2)||(split_cell(i_8_m1(i),j_8(i))==4);
    end   
       
    % Loop to redefine the triangles in the selected quadrants
    element_to_delete=[];
    Cell_to_visit=find(sum(Quadrant_8_no_split,2)>0);
    
    corner_1 =            (i_8r(Cell_to_visit)  )+nb_x*(j_8r(Cell_to_visit)  -1);
    corner_2 =            (i_8r(Cell_to_visit)+1)+nb_x*(j_8r(Cell_to_visit)  -1);
    corner_3 =            (i_8r(Cell_to_visit)+1)+nb_x*(j_8r(Cell_to_visit)+1-1);
    corner_4 =            (i_8r(Cell_to_visit)  )+nb_x*(j_8r(Cell_to_visit)+1-1);
    center   =  nb_x*nb_y+(i_8r(Cell_to_visit)  )+nb_x*(j_8r(Cell_to_visit)  -1);
    mid_down =2*nb_x*nb_y+(i_8r(Cell_to_visit)  )+nb_x*(j_8r(Cell_to_visit)  -1);
    mid_up   =2*nb_x*nb_y+(i_8r(Cell_to_visit)  )+nb_x*(j_8r(Cell_to_visit)+1-1);
    mid_left =3*nb_x*nb_y+(i_8r(Cell_to_visit)  )+nb_x*(j_8r(Cell_to_visit)  -1);
    mid_right=3*nb_x*nb_y+(i_8r(Cell_to_visit)+1)+nb_x*(j_8r(Cell_to_visit)  -1);
    
    local_start=0;
    for i=1:length(Cell_to_visit)
        
        
        new_start=end_8+1;
        tmp_start=start_8+(Cell_to_visit(i)-1)*8;
        for quadrant=1:4,
            
            switch quadrant
                case 1
                    sub_corner=[corner_1(i),mid_down(i),corner_2(i),center(i)];
                case 2
                    sub_corner=[corner_2(i),mid_right(i),corner_3(i),center(i)];
                case 3
                    sub_corner=[corner_3(i),mid_up(i),corner_4(i),center(i)];
                case 4
                    sub_corner=[corner_4(i),mid_left(i),corner_1(i),center(i)];
            end
            
            if(Quadrant_8_no_split(Cell_to_visit(i),quadrant))
                element.num_node(new_start+local_start,1)=sub_corner(1);
                element.num_node(new_start+local_start,2)=sub_corner(3);
                element.num_node(new_start+local_start,3)=sub_corner(4);
                element.depth   (new_start+local_start,1)=depthr(corner_1(i));
                local_start=local_start+1;
            else
                element.num_node(new_start+local_start,1)=sub_corner(1);
                element.num_node(new_start+local_start,2)=sub_corner(2);
                element.num_node(new_start+local_start,3)=sub_corner(4);
                element.depth   (new_start+local_start,1)=depthr(corner_1(i));
                local_start=local_start+1;
                
                element.num_node(new_start+local_start,1)=sub_corner(2);
                element.num_node(new_start+local_start,2)=sub_corner(3);
                element.num_node(new_start+local_start,3)=sub_corner(4);
                element.depth   (new_start+local_start,1)=depthr(corner_1(i));
                local_start=local_start+1;
            end
        end
        element_to_delete=[element_to_delete;(tmp_start+[0:7])'];
    end
    disp('Correct the transition between split_8 and split_4 or split_2')
    disp('Number of elements to delete')
    length(element_to_delete)
    element.num_node(element_to_delete,:)=[];
    element.depth(element_to_delete,:)=[];
end

%---------- Final information --------------
% Final number of elements
mesh.Ne=size(element.num_node,1);

% Information on the nodes
mesh.Nn=4*nb_x*nb_y;
mesh.node.num = 1:mesh.Nn;
mesh.node.lat   = [qlatr(:);platr(:);vlatr(:);ulatr(:)]';
mesh.node.lon   = [qlonr(:);plonr(:);vlonr(:);ulonr(:)]';

disp('loading elements finished')

return

