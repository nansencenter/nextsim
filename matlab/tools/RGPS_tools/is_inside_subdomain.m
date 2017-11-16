function [ inside_subdomain ] = is_inside_subdomain( limit_domain, x_target, y_target )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nb_subdomains=length(limit_domain);
inside_subdomain=zeros(length(x_target),1);

for j=1:nb_subdomains
    sign_cross_product=zeros(length(x_target),1);
    tmp_subdomain=limit_domain(j);
    nb_edges_polygon=length(tmp_subdomain.x)-1;
    for k=1:nb_edges_polygon,
        x1=tmp_subdomain.x(k);
        x2=tmp_subdomain.x(k+1);
        y1=tmp_subdomain.y(k);
        y2=tmp_subdomain.y(k+1);
        sign_cross_product=sign_cross_product+sign((x_target-x1).*(y2-y1)-(y_target-y1).*(x2-x1));
    end
    inside_subdomain=inside_subdomain+(abs(sign_cross_product)==nb_edges_polygon);
end

if(nb_subdomains==0)
    inside_subdomain=inside_subdomain+1;
end

end

