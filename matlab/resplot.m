function resplot(step)

clearvars -except step;

% mc_script=['mu',num2str(step)];
% mx_script=['mx',num2str(step)];
% my_script=['my',num2str(step)];
% 
% run(mc_script);
% run(mx_script);
% run(my_script);
% 
% var_mc=eval(['var_', mc_script]);
% var_mx=eval(['var_', mx_script]);
% var_my=eval(['var_', my_script]);
% 
% [nr,nc]= size(var_mc);
% c=reshape(var_mc,[3,nr/3]);
% x=reshape(var_mx,[3,nr/3]);
% y=reshape(var_my,[3,nr/3]);
% 
% patch(x,y,c,'FaceColor','flat','EdgeColor','none')
% caxis([min(var_mc), max(var_mc)])

field='Velocity';
%field='mld';
%field='Concentration';
%field='Thickness';
%field='Wind';
%field='Ocean';
%field='Vair_factor';
%field='Voce_factor';
%field='Damage';
[mesh_out,data_out] = neXtSIM_bin_revert('', step);

% reshape
var_mx=mesh_out.Nodes_x(mesh_out.Elements);
var_my=mesh_out.Nodes_y(mesh_out.Elements);

[nr,nc]= size(var_mx);
Ne=nr/3;
Nn=length(mesh_out.Nodes_x);
x=reshape(var_mx,[3,Ne]);
y=reshape(var_my,[3,Ne]);

field_tmp=data_out.(field);
length(field_tmp)
if(length(field_tmp)==Ne)
    c{1}=[field_tmp,field_tmp,field_tmp]';
elseif(length(field_tmp)==2*Nn)
    var_mc=field_tmp(mesh_out.Elements);
    c{1}=reshape(var_mc,[3,Ne]);
    var_mc=field_tmp(mesh_out.Elements+Nn);
    c{2}=reshape(var_mc,[3,Ne]);
elseif(length(field_tmp)==Nn)
    var_mc=field_tmp(mesh_out.Elements);
    c{1}=reshape(var_mc,[3,Ne]);
else
    error('Not the right dimensions')
end

for i=1:length(c)
    figure
    patch(x,y,c{i},'EdgeColor','none')
    min(min(c{i}))
    max(max(c{i}))
    caxis([min(min(c{i})), max(max(c{i}))])
    colorbar
end

end
