function sortmst(path2data)
% SORTMST Sorts master records according to the 'From' field
%
% Note: most of the calculations in Analyze rely on sorted mst.

if nargin < 1 || isempty(path2data),  path2data = '.\data\TAQ'; end

% Open matlabpool
if matlabpool('size') == 0
    matlabpool('open', 4, 'AttachedFiles',{'.\utils\poolStartup.m'})
end

dd  = dir(fullfile(path2data,'*.mat'));
N   = numel(dd);
parfor f = 1:N
    disp(f)
    % Load data
    fname = fullfile(path2data,dd(f).name);
    s     = load(fname);
    if ~issorted(s.mst.From)
        disp('Sorting')
        % Sort
        s.mst = sortrows(s.mst,'From');
        % Save
        fileattrib(fname,'+w')
        save(fname, '-struct','s')
        fileattrib(fname,'-w +h')
    end
end
end