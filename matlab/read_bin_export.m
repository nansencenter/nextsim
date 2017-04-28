function out=read_bin_export(file_info,file_data)
%Getting info
% disp(['Read: ' file_info])
[fileID, mssg] = fopen(file_info,'r');
if fileID < 0, error([file_info ': ' mssg]), end

line1 = strsplit(fgetl(fileID));
fseek(fileID,0,'bof');%go back to start of file
N     = length(line1);
fmt   = '%s %s';
if N>=3; fmt   = [fmt,' %d']; end % length of record
if N>=4; fmt   = [fmt,' %f']; end % min
if N>=5; fmt   = [fmt,' %f']; end % max
info_tmp=textscan(fileID,fmt);
fclose(fileID);

data_names=info_tmp{1};
data_types=info_tmp{2};

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

    N_data     = fread(fileID,1,'int');
    position   = ftell(fileID);
    if strcmp('Time',field{i})
        if strcmp(prec,'double')
           Time   = fread(fileID,N_data,prec);
           if (Time>0&Time<1e-12) | Time>1e30
              %%probably read wrong, try again with single
              fseek(fileID,position,'bof');
              Time   = fread(fileID,N_data,'float');
           end
        elseif strcmp(prec,'single')
           Time   = fread(fileID,N_data,prec);
           if (Time>0&Time<1e-12) | Time>1e30
              %%probably read wrong, try again with double
              fseek(fileID,position,'bof');
              Time   = fread(fileID,N_data,'double');
           end
        end
        out.(field{i}) = Time + datenum(1900,1,1);
    else
        out.(field{i}) = fread(fileID,N_data,prec);
    end
end
fclose(fileID);
end
