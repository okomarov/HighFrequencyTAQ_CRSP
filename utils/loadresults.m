function loadresults(name,outname)
if nargin == 1, outname = name; end
resdir  = '.\results';
files   = dir(fullfile(resdir, sprintf('*%s.mat', name)));
[~,idx] = max(files.datenum); % Most recent
s       = load(fullfile(resdir,files(idx).name));
assignin('caller', outname, s.(char(fieldnames(s))));
end