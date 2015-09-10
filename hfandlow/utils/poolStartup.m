function p = poolStartup(varargin)
% POOLSTARTUP Opens pool an autoreferences itself to attach a system call on each worker

% Check if debug mode
pos = find(strcmpi(varargin, 'debug'));
if isempty(pos)
    isdebug = false;
else
    isdebug = varargin{pos+1};
    varargin(pos:pos+1) = [];
end

% Open a pool 
if nargin > 0
    if isempty(gcp('nocreate')) && ~isdebug
        p = parpool(varargin{:});
    end
else
    % Worker priority in win to low
    system(sprintf('wmic process where processid="%d" CALL setpriority 64',feature('getpid')))
end
end