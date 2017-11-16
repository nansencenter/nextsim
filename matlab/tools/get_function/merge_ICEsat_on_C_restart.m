function merge_ICEsat_on_C_restart(saved_simul_in)


simul_in=read_simul_in(saved_simul_in); 
time_init=simul_in.time_init;

%time_init=datenum('15-Oct-2007 12:00:00');%simul_in.time_init;
    
% check if the choosen period is covered by the simulation
dirname='restart';
step=0;
[mesh_out,data_out] = neXtSIM_bin_revert(dirname,[], step);
field_info  = [ dirname '/field_' num2str(step) '.dat'];

%reshape
var_mx=mesh_out.Nodes_x(mesh_out.Elements);
var_my=mesh_out.Nodes_y(mesh_out.Elements);
nr = size(var_mx,1);
Ne=nr/3;
Nn=length(mesh_out.Nodes_x);
x=reshape(var_mx,[3,Ne]);
y=reshape(var_my,[3,Ne]);

element.x=mean(x,1)/1000;
element.y=mean(y,1)/1000;
%time_init=datenum('01-Mar-2008 00:00:00')
time_init=datenum(simul_in.time_init)
h_ICEsat=get_icethickICEsat(time_init,element);

h_init=h_ICEsat;
h_model=data_out.M_thick';
c_model=data_out.M_conc';
%h_model=data_out.Thickness';
%c_model=data_out.Concentration';

h_init=h_model;
hmin_confidence=2;
mask=((~isnan(h_ICEsat)).*(h_model>0))>0;
confidence_ICEsat=min(1.,h_ICEsat(mask)./hmin_confidence.*c_model(mask));
h_init(mask)=confidence_ICEsat.*h_ICEsat(mask)+(1-confidence_ICEsat).*h_model(mask);

% reshaping
h_ICEsat=[h_ICEsat;h_ICEsat;h_ICEsat];
h_model=[h_model;h_model;h_model];
h_new=[h_init;h_init;h_init];


%ice mask and water mask extraction
mask=h_init;
mask_ice=find(mask>0);
mask_water=find(mask==0);

figure
patch(x/1000,y/1000,h_new,'FaceColor','flat','EdgeColor','none')
caxis([0 4])
figure
patch(x/1000,y/1000,h_model,'FaceColor','flat','EdgeColor','none')
caxis([0 4])
figure
patch(x/1000,y/1000,h_ICEsat,'FaceColor','flat','EdgeColor','none')
caxis([0 4])

data_out.M_thick=h_init;
data_out.M_damage(data_out.M_damage==0.)=0.5;

write_bin_export(data_out,field_info,'field_0_modified.bin');

end
