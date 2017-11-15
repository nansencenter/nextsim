function failed = resplot(field,step,dirname,plot_options)
%% CALL: failed = resplot(field,step,dirname,plot_options)

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

% plot_options = struct eg
%     visible        : 1
%   plot_grid        : 1
%  vector_components : []

failed   = 0;
if ~exist('dirname','var'), dirname='.'; end
if ~exist('step','var'); step='init'; end

%other plot options:
if exist('plot_options','var')
   if ~isstruct(plot_options)
      % backwards compatibility
      plot_grid   = plot_options;
   else
      flds  = fieldnames(plot_options);
      for j=1:length(flds)
         fld   = flds{j};
         cmd   = [fld,' = plot_options.',fld,';'];
         eval(cmd);
      end
   end
end

if ~exist('plot_grid','var'),          plot_grid = 0; end;
   % If not zero the mesh lines are ploted. If zoomed out only the mesh lines will be visible
if ~exist('visible','var'),            visible = 1; end;
   % we display the figure on the screen (may be set to 0 when generating a large amount of figures)
if ~exist('vector_components','var'),  vector_components = []; end;
   % if a vector: vector_components is [] (plot all components) or vector of components to plot,
   % where i=1 is x-component, i=2 is y-component, i=3 is modulus
if ~exist('no_error','var'), no_error = 0; end;


[mesh_out,data_out] = neXtSIM_bin_revert(dirname,[], step);

% reshape
var_mx=mesh_out.Nodes_x(mesh_out.Elements);
var_my=mesh_out.Nodes_y(mesh_out.Elements);

if(isfield(data_out,'M_dirichlet_flags'))
    
    dirichlet_x=mesh_out.Nodes_x(data_out.M_dirichlet_flags+1);
    dirichlet_y=mesh_out.Nodes_y(data_out.M_dirichlet_flags+1);
    %figure
    %plot(dirichlet_x,dirichlet_y,'.')
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
        if ~quiet
           disp(['Available fields are: ' fieldnames(data_out)']);
        end
        if no_error
           failed = 1;
           return;
        else
           rethrow(err);
        end
    end
else
    field_tmp=zeros(Ne,1);
    plot_grid=true;
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


if visible
   visible  = 'on';
else
   visible  = 'off';
end

if length(c)>1
   %vector
   if ~isempty(vector_components)
      c  = c(vector_components);
   end
end
for i=1:length(c)
    figure('visible',visible);
    if(plot_grid)
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
