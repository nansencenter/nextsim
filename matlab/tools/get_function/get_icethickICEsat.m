function h_init=get_icethickICEsat(time_init,element)

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
    if((Period_start{i}<=time_init)&&(Period_end{i}>=time_init))
        selected_campaign=i;
        ICESAT_campaign=Campaign{i};
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
    
    H1=IceSat.Th/100;
    x=IceSat.x;
    y=IceSat.y;

    %To change the no ice mask to zero 
    H1(H1==-1/100)=NaN;

    % Now the interpolation onto the nodes of our finite element mesh
    f=(isnan(H1)==0);
    h_init=griddata(x(f),y(f),H1(f),element.x,element.y,'linear');
    
    % Put NaN to points South of Fram Strait 
    x1=-1500;
    y1=-1500;
    x2=1000;
    y2=-550;
    South_Fram_Strait=(((element.y-y1)./(element.x-x1))<((y2-y1)./(x2-x1))).*((element.x-x1)>0);
    h_init(South_Fram_Strait==1)=NaN;
    
%     % Use the nearest neigbour with a decreasing factor depending on the
%     % distqnance to observation
%     f_interp=(isnan(h_init)==1);
%     h_init(f_interp)=griddata(x(f),y(f),H1(f),element.x(f_interp),element.y(f_interp),'nearest');
%     x_data=element.x;
%     y_data=element.y;
%     x_data(f_interp)=griddata(x(f),y(f),x(f),element.x(f_interp),element.y(f_interp),'nearest');
%     y_data(f_interp)=griddata(x(f),y(f),y(f),element.x(f_interp),element.y(f_interp),'nearest');
%     distance_to_data=hypot(x_data-element.x,y_data-element.y);
%     half_radius=200; %km
%     h_init=h_init./(1+distance_to_data/half_radius);

else
    error(['no ICESAT campaign corresponds are available for that date:' datestr(time_init)])
end