function simul_in=read_simul_in(saved_simul_in,DO_DISP)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Example: simul_in=read_simul_in('log_simul')

if ~exist('DO_DISP','var')
   DO_DISP  = 1;
end

fields   = {'time_init',...
            'duration',...
            'output_per_day',...
            'undamaged_time_relaxation_sigma',...
            'exponent_relaxation_sigma',...
            'young',...
            'mesh_filename'};

Nfld  = length(fields);
for j=1:Nfld
   fld   = fields{j};
   lin   = get_value(saved_simul_in,['simul.',fld],DO_DISP);
   nlin  = str2num(lin);
   if ~isempty(nlin)
      % see if it's a number
      simul_in.(fld) = nlin;
   else
      % leave as a string
      simul_in.(fld) = lin;
   end
end

try
    simul_in.time_init=datenum(simul_in.time_init,'yyyy-mm-ddHH:MM:SS');
catch
    try
        simul_in.time_init=datenum(simul_in.time_init,'yyyy-mm-dd');
    catch err
        warning('Couldn''t read time from nextsim.log - check the file!')
        simul_in.time_init=datenum(1900,1,1);
    end
end

end

%-----------------------
function A=get_value(saved_simul_in,look_for,DO_DISP)

fid = fopen(saved_simul_in);
tline = fgetl(fid);
while ischar(tline)
    
    if(length(tline)>=length(look_for))
       if(strcmp(tline(1:length(look_for)),look_for))
          if DO_DISP
             disp(tline);
          end
       
          A = sscanf(tline(length(look_for)+1:end), '%s');
       end
    end
    tline = fgetl(fid);
end

fclose(fid);

end

