function [mesh_out,data_out] = neXtSIM_bin_revert(dirname, proc, step,SHOW_INFO)%,varargin)

%Sylvain: Created 20160105 to take a binary file and revert it to a slim_out.
%Only works for binary files created by neXtSIM
%example,
%>> neXtSIM_bin_revert('',[],0)
%>> neXtSIM_bin_revert('',[],'0')
%>> neXtSIM_bin_revert('',[],'init')
%>> neXtSIM_bin_revert('',[],'final')

%defaults:
% save_flag = 0; %don't save
% 
% nVarargs = length(varargin);
% if nVarargs >= 1, save_flag = varargin{1}; end
% if nVarargs >= 2, error('Too many inputs'), end

if ~exist('dirname','var'); dirname  = './'; end
if ~exist('proc','var'); proc  = []; end
if ~exist('step','var'); step  = 'final'; end
if ~exist('SHOW_INFO','var'); SHOW_INFO = 0; end

if(~isempty(dirname)&& dirname(end)~='/')
    dirname=[dirname, '/'];
end

if(isempty(step));
   step = 'final';
elseif ~ischar(step)
   step  = num2str(step);
end

if(isempty(proc))
    field_info  = [ dirname 'field_' step '.dat'];
    field_data  = [ dirname 'field_' step '.bin'];
    mesh_info   = [ dirname 'mesh_'  step '.dat'];
    mesh_data   = [ dirname 'mesh_'  step '.bin'];
else
    if ischar(proc)
       sproc   = proc;
    else
       sproc   = num2str(proc);
    end
    field_info  = [ dirname 'field_' sproc '_' step '.dat'];
    field_data  = [ dirname 'field_' sproc '_' step '.bin'];
    mesh_info   = [ dirname 'mesh_'  sproc '_' step '.dat'];
    mesh_data   = [ dirname 'mesh_'  sproc '_' step '.bin'];
end


%Getting mesh
mesh_out=read_bin_export(mesh_info,mesh_data,SHOW_INFO);

%Getting field
if(exist(field_info))
    data_out=read_bin_export(field_info,field_data,SHOW_INFO);
else
    data_out=[];
end

return
