function resplot(step)

clearvars -except step;

mc_script=['mu',num2str(step)];
mx_script=['mx',num2str(step)];
my_script=['my',num2str(step)];

run(mc_script);
run(mx_script);
run(my_script);

var_mc=eval(['var_', mc_script]);
var_mx=eval(['var_', mx_script]);
var_my=eval(['var_', my_script]);

[nr,nc]= size(var_mc);
c=reshape(var_mc,[3,nr/3]);
x=reshape(var_mx,[3,nr/3]);
y=reshape(var_my,[3,nr/3]);

patch(x,y,c,'FaceColor','flat','EdgeColor','none')
caxis([min(var_mc), max(var_mc)])
end
