function out=read_bin_export(file_info,file_data)
%Getting info
% disp(['Read: ' file_info])
[fileID, mssg] = fopen(file_info,'r');
if fileID < 0, error([file_info ': ' mssg]), end
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