function sortmst(path2data)
% SORTMST Sorts master records according to the 'From' field
%
% Note: most of the calculations in Analyze rely on sorted mst.

if nargin < 1 || isempty(path2data),  path2data = '.\data\TAQ'; end

setupemail

% Open matlabpool
if isempty(gcp('nocreate'))
    parpool(4, 'AttachedFiles',{'.\utils\poolStartup.m'})
end

files = dir(fullfile(path2data,'*.mat'));
N     = numel(files);
tic
parfor f = 1:N
    disp(f)
    % Load data
    fname = fullfile(path2data,files(f).name);
    s     = load(fname);
    if ~issorted(s.mst.From)
        disp('Sorting')
        % Sort
        s.mst = sortrows(s.mst,'From');
        % Save
        fileattrib(fname,'+w')
        mysave(fname, s)
    end
    fileattrib(fname,'-w +h')
end
% Cleanup and feedback
delete(gcp)
message = sprintf('Task ''sortmst'' terminated in %s',sec2time(toc));
disp(message)
sendmail('o.komarov11@imperial.ac.uk', message,'');
rmpref('Internet','SMTP_Password')
end

function mysave(fname, s)
save(fname, '-struct','s')
end