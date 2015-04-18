function out = loadresults(name)
% LOADRESULTS Load results from .\results
%
%   LOADRESULTS(NAME) Loads the latest '*NAME.mat'

resdir  = '.\results';
if isempty(strfind(name,'.mat'))
    name = [name, '.mat'];
end
files   = dir(fullfile(resdir, ['*', name]));
if isempty(files)
    error('loadresults:nofile','No files matching ''%s'' found.',name)
end
[~,idx] = max([files.datenum]); % Most recent
s       = load(fullfile(resdir,files(idx).name));
out     = s.(char(fieldnames(s)));
if isa(out, 'dataset'), out = dataset2table(out); end
end