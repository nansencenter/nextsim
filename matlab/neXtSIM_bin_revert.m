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

if(isempty(proc))
    field_info  = [ dirname '/field_' num2str(step) '.dat'];
    field_data  = [ dirname '/field_' num2str(step) '.bin'];
    mesh_info   = [ dirname '/mesh_' num2str(step) '.dat'];
    mesh_data   = [ dirname '/mesh_' num2str(step) '.bin'];
else
    field_info  = [ dirname '/field_' num2str(proc) '_' num2str(step) '.dat'];
    field_data  = [ dirname '/field_' num2str(proc) '_' num2str(step) '.bin'];
    mesh_info   = [ dirname '/mesh_' num2str(proc) '_' num2str(step) '.dat'];
    mesh_data   = [ dirname '/mesh_' num2str(proc) '_' num2str(step) '.bin'];
end


%Getting mesh
mesh_out=read_bin_export(mesh_info,mesh_data);

%Getting field
data_out=read_bin_export(field_info,field_data);

end

function out=read_bin_export(file_info,file_data)
%Getting info
fileID = fopen(file_info,'r');
info_tmp=textscan(fileID,'%s %s');
data_names=info_tmp{1};
data_types=info_tmp{2};
fclose(fileID);

field = cell(numel(data_names));
fileID = fopen(file_data,'r');
for i = 1:numel(data_names)
    field{i} = data_names{i};
    %if data_types(i)==8
    %   prec = 'double';
    % else
    %   prec = 'int';
    % end

    if strcmp(data_types(i),'double') || strcmp(data_types(i),'8') % check against '8' for backwards compatiblity
	  prec = 'double';
	elseif strcmp(data_types(i),'float')
	  prec = 'single';
	else
	  prec = 'int';
    end

    N_data=fread(fileID,1,'int');
    if strcmp('Time',field{i})
        out.(field{i}) = fread(fileID,N_data,prec) + datenum(1900,1,1);
    else
        out.(field{i}) = fread(fileID,N_data,prec);
    end
end
fclose(fileID);
end
