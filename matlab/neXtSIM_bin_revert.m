function [mesh_out,data_out] = neXtSIM_bin_revert(dirname, proc, step)%,varargin)

%Sylvain: Created 20160105 to take a binary file and revert it to a slim_out.
%Only works for binary files created by neXtSIM
%example,
%>> neXtSIM_bin_revert('',0)

%defaults:
% save_flag = 0; %don't save
% 
% nVarargs = length(varargin);
% if nVarargs >= 1, save_flag = varargin{1}; end
% if nVarargs >= 2, error('Too many inputs'), end

if(~isempty(dirname)&& dirname(end)~='/')
    dirname=[dirname, '/'];
end

if(isempty(proc))
    field_info  = [ dirname 'field_' num2str(step) '.dat'];
    field_data  = [ dirname 'field_' num2str(step) '.bin'];
    mesh_info   = [ dirname 'mesh_' num2str(step) '.dat'];
    mesh_data   = [ dirname 'mesh_' num2str(step) '.bin'];
else
    field_info  = [ dirname 'field_' num2str(proc) '_' num2str(step) '.dat'];
    field_data  = [ dirname 'field_' num2str(proc) '_' num2str(step) '.bin'];
    mesh_info   = [ dirname 'mesh_' num2str(proc) '_' num2str(step) '.dat'];
    mesh_data   = [ dirname 'mesh_' num2str(proc) '_' num2str(step) '.bin'];
end


%Getting mesh
mesh_out=read_bin_export(mesh_info,mesh_data);

%Getting field
if(exist(field_info))
    data_out=read_bin_export(field_info,field_data);
else
    data_out=[];
end

end

