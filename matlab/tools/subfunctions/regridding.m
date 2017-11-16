function [ebox]=regridding(order,xmin,ymin,min_box_size,cells_x,cells_y,cells_stream,cells_date,cells_deltat,date)

%function []=regridding()
%Auteur: Sylvain Bouillon, June 2013
%
%GOAL: regridding to ease masking and coarse-graining
%
%OUTPUT:

ebox.mask=[];
ebox.no_overlap=[];
ebox.full=[];

if(isempty(cells_x) || isempty(order))
    return;
end

if(nargin<8)
    cells_stream=[];
    cells_date=[];
    cells_deltat=[];
    date=0;
    deltat=0;
end

ratio_dom_resol = 2^order;

ebox.xmin=xmin;
ebox.ymin=ymin;
ebox.min_box_size=min_box_size;
ebox.ratio_dom_resol=ratio_dom_resol;

% Initialisation
% RGPS cells can be superposed as coming from overlaping stream
% We select only RGPS cells coming from a unique stream for each box at the lowest resol 
ebox.mask=zeros(ratio_dom_resol,ratio_dom_resol);
ebox.full=cell(ratio_dom_resol,ratio_dom_resol);
ebox.no_overlap=cell(ratio_dom_resol,ratio_dom_resol);

indices_x=floor((cells_x-xmin)/min_box_size)+1;
indices_y=floor((cells_y-ymin)/min_box_size)+1;

indices_tot=(indices_x-1)*ratio_dom_resol+indices_y;

% define a list of unique box indices
[indice_unique,i_unique_tot,i_tot_unique] = unique(indices_tot);

% We count the occurence of each box indices to speed up the search
if length(indice_unique) == 1, counts = length(indices_tot);
else counts = hist(indices_tot,indice_unique); end

% We sort the list of cell indices to speed up the search
[ind_sorted,sort_ind]=sort(i_tot_unique);
indices_list=1:length(indices_tot);
indices_list_sorted=indices_list(sort_ind);

counter_tmp=1;
for i=1:length(indice_unique)
    % old very expensive way to search the indices
    % ind=find(indices_tot==indice_unique(i)); 
    
    % new pretty and efficient way to search the indices thanks to unique()
    new_counter_tmp=counter_tmp+counts(i);
    ind=indices_list_sorted(counter_tmp:new_counter_tmp-1);
    counter_tmp=new_counter_tmp;
    
    ind_local_x=floor(indice_unique(i)/ratio_dom_resol)+1;
    ind_local_y=indice_unique(i)-(ind_local_x-1)*ratio_dom_resol;
    
    if(ind_local_x>0 && ind_local_y>0 && ind_local_x<ratio_dom_resol+1 && ind_local_y<ratio_dom_resol+1)
        
        % We save this information for the comparison
        nb_elements=length(ind);
        ebox.mask(ind_local_x,ind_local_y) = 1;
        ebox.full{ind_local_x,ind_local_y} = [ebox.full{ind_local_x,ind_local_y},ind];
                
        % We now select only the data coming from the same date and stream
        % allows only data coming from the same date
        if(~isempty(cells_date))
            [date_dist,date_ind] = min(abs(cells_date(ind)-date));
            date_selected = cells_date(ind(date_ind));
            deltat_selected = cells_deltat(ind(date_ind));
            ind= ind(find(((cells_date(ind)+cells_deltat(ind)/2)==(date_selected+deltat_selected/2))));
        end
        % allows only data coming from the same stream to avoid surperposed RGPS streams
        if(~isempty(cells_stream) && length(ind)>1)
            stream_selected = cells_stream(ind(1));
            ind= ind(cells_stream(ind)==stream_selected);
        end
        
        if(isempty(ind))
            error('isempty_ind')
        end
        
        nb_elements=length(ind);
        ebox.no_overlap{ind_local_x,ind_local_y} = [ebox.no_overlap{ind_local_x,ind_local_y},ind];
    end
end
