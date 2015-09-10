function out = loadresults(name, resdir)
% LOADRESULTS Load results from .\results
%
%   LOADRESULTS(NAME) Loads the latest '*NAME.mat'

if nargin < 2
    resdir  = '.\results';
end

if isempty(strfind(name,'.mat'))
    name = [name, '.mat'];
end
% Most recent
files   = dir(fullfile(resdir, ['*', name]));
if isempty(files)
    error('loadresults:nofile','No files matching ''%s'' found.',name)
end
[~,idx] = max([files.datenum]); 
s       = load(fullfile(resdir,files(idx).name));
% Load whole structure if multiple fields
fnames  = fieldnames(s);
if numel(fnames) == 1
    out = s.(fnames{1});
else
    out = s;
end
% Convert datasets to table
if isa(out, 'dataset'), out = dataset2table(out); end
end