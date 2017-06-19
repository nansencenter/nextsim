function simul_in=read_simul_in(saved_simul_in,DO_DISP)

if ~exist('saved_simul_in','var'); saved_simul_in  = 'nextsim.log'; end
if ~exist('DO_DISP','var'); DO_DISP  = 0; end

%treat this differently, since can have multiple config files
simul_in.info.config_files   = {};

fid   = fopen(saved_simul_in);

%%read in rest of variables:
while ~feof(fid)
   [x,name] = read_next(fid);
   if ~isempty(name)
      if strcmp(name,'config_files')
         simul_in.info.config_files{end+1}  = x;
      else
         eval(['simul_in.',name,' = x;']);
      end
   end
end
fclose(fid);

if DO_DISP
   fields   = fieldnames(simul_in);
   for n=1:length(fields)
      fld   = fields{n};
      disp([fld,':']);
      disp(simul_in.(fld));
   end
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,name]  = read_next(fid)
%% read next line in text file

lin   = strtrim(fgets(fid));  %% trim leading white space
lin2  = strsplit(lin);        %%split using spaces
name  = lin2{1};              %%get 1st thing in line (as string)

% blank line
if strcmp(name,'')
   x     = [];
   name  = [];
   return
end

% comment
if strcmp(lin(1),'#')
   x     = [];
   name  = [];
   return
end


% ===============================================================
% special cases
% - config files
if~isempty(strfind(lin,']='))
   name  = 'config_files';
   x     = strsplit(lin,']=');
   x     = x{2};
   return
end

% - spaces in names)
if length(lin2)>2
   if strcmp('Git',lin2{1}) & strcmp('revision',lin2{2})
      name  = 'info.Git_revision';
      x     = strtrim(lin(length(name)-4:end));
      return
   end
   if strcmp('Build',lin2{1}) & strcmp('date',lin2{2})
      name  = 'info.Build_date';
      x     = strtrim(lin(length(name)-4:end));
      return
   end
end

% - name too long
if length(name)==length(lin);
   x     = NaN;
   msg   = ['error in reading simul_in: ',name ' = NaN'];
   warning(msg);
   return
end

% names that are dates (convert to datenum)
time_fields   = {'simul.time_init'};%list of names that are dates
time_field    = [];
for n=1:length(time_fields)
   if strcmp(time_fields{n},name)
      try
         x  = datenum(strtrim(lin(length(name)+1:end)));
      catch
         warning(['Couldn''t read time (',name,') from nextsim.log - check the file!'])
         x  = datenum(1900,1,1);
      end
      return;
   end
end
% ===============================================================

% proper variable
xs          = strtrim(lin(length(name)+1:end));
[x,status]  = str2num(xs);
if status==0
   %% numerical conversion failed
   %% - leave as string
   x  = xs;
end
if isempty(strfind(name,'.'))
   name  = ['info.',name];
end

% remove '-' from name (can't be a field name for a structure)
name  = strrep(name,'-','_');

% rename some fields to make them more informative:
if strcmp(name,'info.C++')
   % also can't have '+' in name (can't be a field name for a structure)
   name  = 'info.Cpp_compiler';
end
if strcmp(name,'info.C')
   name  = 'info.C_compiler';
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
