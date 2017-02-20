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
p.CaseSensitive = false;
addParameter(p,'Debug',false);
addParameter(p,'NumCores',4);
addParameter(p,'CacheByFile',cell(N,1));
addParameter(p,'InputToFun',{});
addParameter(p,'DisplayProgress',true);
p.parse(varargin{:});

CacheByFile  = p.Results.CacheByFile;
InputToFun   = p.Results.InputToFun;
DispProgress = p.Results.DisplayProgress;

% Setup
poolStartup(p.Results.NumCores, p.Results.Debug)
if ~p.Results.Debug
    setupemail();
end

try
    tic
    parfor f = 1:N
        disp(f)
        fname  = fullfile(path2data, dd(f).name);
        s      = taq.util.loadFull(fname);
        inputs = {CacheByFile{f,:}, f, InputToFun{:}};
        res{f} = fun(s, inputs{:});
    end
    res = cat(1,res{:});

    % Export results and notify
    filename = sprintf('%s_%s.mat',outname,datestr(now,'yyyymmdd_HHMM'));
    save(filename, 'res')
    stack    = dbstack();
    message  = sprintf('Task ''%s'' terminated in %s',stack(end).name,sec2time(toc));
    disp(message);
    if ~p.Results.Debug
        sendmail('o.komarov11@imperial.ac.uk', message,'');
    end
catch err
    filename = fullfile(fileparts(outname), sprintf('err_%s.mat',datestr(now,'yyyymmdd_HHMM')));
    save(filename,'err')
    if ~p.Results.Debug
        stack = dbstack();
        sendmail('o.komarov11@imperial.ac.uk',sprintf('ERROR in task ''%s''',stack(end).name), err.message, {filename})
    end
    rethrow(err)
end
end
