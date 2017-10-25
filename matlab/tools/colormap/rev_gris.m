function map = rev_gris(n)
%rev_gris
%
%n=256;
if nargin < 1, n = size(get(gcf,'colormap'),1); end
gris=hot(n);
rev_gris((n:-1:1),:)=gris(-(-n:-1),:);
tmp=rev_gris(:,3); rev_gris(:,3)=rev_gris(:,1); rev_gris(:,1)=tmp;
rev_gris(end,1)=0.98; rev_gris(end,2)=1.0; rev_gris(end,3)=1.0;
map=rev_gris;
