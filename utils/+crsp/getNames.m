function names = getNames(freq)
% GETNAMES Loads the CRSP stock event - name history dataset
persistent names_d names_m

if nargin < 1 || isempty(freq)
    freq = 'monthly';
end

freq = validateFreq(freq);

folder = fullfile(fileparts(mfilename('fullpath')),'data');

switch freq
    case 'monthly'
        if isempty(names_m)
            names_m = loadlatest('msenames',folder);
        end
        names = names_m;
    case 'daily'
        if isempty(names_d)
            names_d = loadlatest('dsenames',folder);
        end
        names = names_d;
    otherwise
        error('Unrecognized FREQ: "%s".',freq)
end
end

function freq = validateFreq(freq)
if ~isrowchar(freq)
    error('crsp:getNames:invalidFreq','FREQ should be a string.')
else
    validstr = {'monthly','daily'};
    pos      = strncmpi(freq, validstr,numel(freq));
    freq     = validstr{pos};
end
end
