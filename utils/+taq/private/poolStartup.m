function p = poolStartup(numcores,isdebug)
% POOLSTARTUP Opens pool an autoreferences itself to attach a system call on each worker

% Open a pool 
if nargin > 0
    if isempty(gcp('nocreate')) && ~isdebug
        p = parpool(numcores, 'AttachedFiles', {'poolStartup.m'});
    end
    
% Worker priority in win to low
else
    system(sprintf('wmic process where processid="%d" CALL setpriority 64',feature('getpid')))
end
end
