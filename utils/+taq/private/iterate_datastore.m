function [res, filename] = iterate_datastore(path2data, outname, fun, varargin)
% ITERATE_DATASTORE Executes specified fun in parallel on the whole filestore (.mat  and .idx files)
%

dd = dir(fullfile(path2data,'*.mat'));
if isempty(dd)
    error('No data found in "%s".', path2data)
end
N   = numel(dd);
res = deal(cell(N,1));

% Parse inputs
p = inputParser();
addParameter(p,'debug',false);
addParameter(p,'numcores',4);
addParameter(p,'cachebyfile',cell(N,1));
addParameter(p,'input2fun',{});
p.parse(varargin{:});

cachebyfile = p.Results.cachebyfile;
input2fun   = p.Results.input2fun;

% Setup
poolStartup(p.Results.numcores, p.Results.debug)
if ~p.Results.debug
    setupemail();
end

try
    tic
    parfor f = 1:N
        disp(f)
        fname  = fullfile(path2data, dd(f).name);
        s      = taq.util.loadFull(fname);
        inputs = {cachebyfile{f,:}, f, input2fun{:}};
        res{f} = fun(s, inputs{:});
    end
    res = cat(1,res{:});

    % Export results and notify
    filename = sprintf('%s_%s.mat',outname,datestr(now,'yyyymmdd_HHMM'));
    save(filename, 'res')
    stack = dbstack();
    message  = sprintf('Task ''%s'' terminated in %s',stack(end).name,sec2time(toc));
    disp(message);
    if ~p.Results.debug
        sendmail('o.komarov11@imperial.ac.uk', message,'');
    end
catch err
    filename = fullfile(fileparts(outname), sprintf('err_%s.mat',datestr(now,'yyyymmdd_HHMM')));
    save(filename,'err')
    if ~p.Results.debug
        stack = dbstack();
        sendmail('o.komarov11@imperial.ac.uk',sprintf('ERROR in task ''%s''',stack(end).name), err.message, {filename})
    end
    rethrow(err)
end
end
