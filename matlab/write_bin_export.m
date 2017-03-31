function write_bin_export(out,file_info,file_bis_data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output precission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int    = 'int32';
double = 'float64';

%Getting info
disp(['Read: ' file_info])
fileID = fopen(file_info,'r');
info_tmp=textscan(fileID,'%s %s');
data_names=info_tmp{1};
data_types=info_tmp{2};
fclose(fileID);

%field = cell(numel(data_names));
fileID_bis = fopen(file_bis_data,'w')

field = data_names;
for i = 1:numel(data_names)   
    %if data_types(i)==8
    %   prec = 'double';
    % else
    %   prec = 'int';
    % end
    if strcmp(data_types(i),'double') || strcmp(data_types(i),'8') % check against '8' for backwards compatiblity
	  prec = 'float64';
	elseif strcmp(data_types(i),'float')
	  %prec = 'single';
      prec = 'float32';
    else
	  prec = 'int32';
    end
    
    N_data=length(out.(field{i}));

   COUNT=fwrite(fileID_bis,N_data,int);
   if(COUNT~=1)
       error('Errtor, not the right precision maybee: ', num2str(COUNT), ' numbers instead of ',num2str(1))
   end
   
   if strcmp('Time',field{i})
       size(out.(field{i}) - datenum(1900,1,1))
       COUNT=fwrite(fileID_bis,out.(field{i}) - datenum(1900,1,1),prec);
   else
       COUNT=fwrite(fileID_bis,out.(field{i}),prec);
   end
   if(COUNT~=N_data)
       error('Errtor, not the right precision maybee: ', num2str(COUNT), ' numbers instead of ',num2str(N_data))
   end
end
fclose(fileID_bis);
end