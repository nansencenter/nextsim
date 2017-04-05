function simul_in=read_simul_in(saved_simul_in)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Example: simul_in=read_simul_in('log_simul')

try
    simul_in.time_init=datenum(get_value(saved_simul_in,'simul.time_init'),'yyyy-mm-ddHH:MM:SS');
catch
    simul_in.time_init=datenum(get_value(saved_simul_in,'simul.time_init'),'yyyy-mm-dd');
end
simul_in.duration=str2num(get_value(saved_simul_in,'simul.duration'));
simul_in.output_per_day=str2num(get_value(saved_simul_in,'simul.output_per_day'));
simul_in.mesh_filename=get_value(saved_simul_in,'simul.mesh_filename');

end

%-----------------------
function A=get_value(saved_simul_in,look_for)

fid = fopen(saved_simul_in);
tline = fgetl(fid);
while ischar(tline)
    
    if(length(tline)>=length(look_for))
        if(strcmp(tline(1:length(look_for)),look_for))
            disp(tline);
        
            A = sscanf(tline(length(look_for)+1:end), '%s');
        end
    end
    tline = fgetl(fid);
end

fclose(fid);

end

