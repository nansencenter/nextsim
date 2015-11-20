clear;
mc;
mx;
my;
[nr,nc]= size(var_mc);
c=reshape(var_mc,[3,nr/3]);
x=reshape(var_mx,[3,nr/3]);
y=reshape(var_my,[3,nr/3]);

patch(x,y,c,'FaceColor','flat','EdgeColor','none')
