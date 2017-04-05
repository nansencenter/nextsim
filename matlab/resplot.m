function resplot(field,step,dir,plot_edge)
%% CALL: resplot(field,step,dir)
if nargin==1, step='init'; end
if nargin<=2, dir=''; end
if nargin<=3, plot_edge=false; end

% clearvars -except step;

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

%field='Velocity';
%field='mld';
%field='Concentration';
%field='Thickness';
%field='SSS';
%field='SST';
%field='Wind';
%field='Ocean';
%field='Vair_factor';
%field='Voce_factor';
%field='Damage';
%field='bathy';
[mesh_out,data_out] = neXtSIM_bin_revert(dir,[], step);

% reshape
var_mx=mesh_out.Nodes_x(mesh_out.Elements);
var_my=mesh_out.Nodes_y(mesh_out.Elements);

if(isfield(data_out,'M_dirichlet_flags'))
    
    dirichlet_x=mesh_out.Nodes_x(data_out.M_dirichlet_flags+1);
    dirichlet_y=mesh_out.Nodes_y(data_out.M_dirichlet_flags+1);
    figure
    plot(dirichlet_x,dirichlet_y,'.')
end

[nr,nc]= size(var_mx);
Ne=nr/3;
Nn=length(mesh_out.Nodes_x);
x=reshape(var_mx,[3,Ne]);
y=reshape(var_my,[3,Ne]);

if(~isempty(field))
    try
        field_tmp=data_out.(field);
    catch err
        disp(['Available fields are: ' fieldnames(data_out)'])
        rethrow(err)
    end
else
    field_tmp=zeros(Ne,1);
    plot_edge=true;
end

if(length(field_tmp)==Ne)
    c{1}=[field_tmp,field_tmp,field_tmp]';
elseif(length(field_tmp)==2*Nn)
    var_mc=field_tmp(mesh_out.Elements);
    c{1}=reshape(var_mc,[3,Ne]);
    var_mc=field_tmp(mesh_out.Elements+Nn);
    c{2}=reshape(var_mc,[3,Ne]);
    c{3}=hypot(c{1},c{2})
elseif(length(field_tmp)==Nn)
    var_mc=field_tmp(mesh_out.Elements);
    c{1}=reshape(var_mc,[3,Ne]);
elseif(length(field_tmp)==1)
    disp([field ' = ' num2str(field_tmp)])
    error(' Not a field ')
else
    error('Not the right dimensions')
end

for i=1:length(c)
    figure
    if(plot_edge)
        patch(x,y,c{i});
    else
        patch(x,y,c{i},'EdgeColor','none')
    min_value=min(min(c{i}))
    max_value=max(max(c{i}))
    if(min_value<max_value)
        caxis([min_value, max_value])
    end
    colorbar
    
    try
        font_size=12;
        textstring=datestr(data_out.Time);
        text(0.55, 0.95,textstring,'units','normalized','BackgroundColor','white','FontSize',font_size,'EdgeColor','k')
    catch e
        display(e)
    end
    set(gca,'DataAspectRatio',[1 1 1], 'Color', [.7 .7 .7])
    end
    
if(~isempty(field))
    % % useful to highligth zero value    
    hold on
    %plot(var_mx(c{1}==0),var_my(c{1}==0),'or')
end
%hold on
%col_fail_id=592575;
%plot(var_mx(col_fail_id+1),var_my(col_fail_id+1),'o')
%axis([var_mx(col_fail_id+1)-1e5,var_mx(col_fail_id+1)+1e5,var_my(col_fail_id+1)-1e5,var_my(col_fail_id+1)+1e5])

end
