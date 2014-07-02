function loadresults(name,outname)
% LOADRESULTS Load results from .\results
%
%   LOADRESULTS(NAME, OUTNAME) Loads the '*NAME.mat' file and assigns to
%                              OUTNAME in the caller ws. If not specified,
%                              OUTNAME = NAME.

if nargin == 1, outname = name; end
resdir  = '.\results';
files   = dir(fullfile(resdir, sprintf('*%s.mat', name)));
[~,idx] = max([files.datenum]); % Most recent
s       = load(fullfile(resdir,files(idx).name));
assignin('caller', outname, s.(char(fieldnames(s))));
end