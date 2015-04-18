function spy = getSpy(freq)
% spy = getSpy(freq)
%
% Note: SPY has permno 84398


% Sampling params
if nargin < 1 || isempty(freq), freq = 5; end
name = sprintf('spy%dm.mat',freq);

% Load
try
    spy = loadresults(name);
% Sample    
catch
    fprintf('%s: sampling SPY at %d min.\n', mfilename, freq)
    opt.grid = (9.5/24:freq/(60*24):16/24)';
    writeto  = '.\results';
    
    % Sample at x min
    [spy, filename] = Analyze('sampleSpy',[],[],[],[],opt);
    
    % Rename to append the sampling frequency
    newfullname = fullfile(writeto, name);
    movefile(fullfile(writeto,filename), newfullname);
end
end