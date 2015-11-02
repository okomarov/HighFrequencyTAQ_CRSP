function [res, filename] = blockprocess(fhandle, projectpath, varnames, cached, path2data, debug, poolcores, varargin)

% BLOCKPROCESS Executes specified fun in parallel on the whole database (all .mat files)
%
%   BLOCKPROCESS(FUN, VARNAMES) FUN should a string with the name of one of
%                          the following sub-functions:
%                               - 'dailystats'
%                               - 'badprices'
%                               - 'avgtimestep'
%                          VARNAMES is a cell-array of strings (or string)
%                          with the VarNames of the dataset with the results
%
%   BLOCKPROCESS(..., PATH2DATA) If you wanna use other than '.\data\TAQ\T*.mat'
%                           files (default), then specify a different
%                           PATH2DATA with the initial pattern of the name,
%                           e.g '.\data\TAQ\sampled\S5m_*.mat'
%
%   BLOCKPROCESS(..., CACHED) Some FUN might require pre-cached results which where
%                        run on the whole database.
%                        Check the specific sub-function for the format of the
%                        needed CACHED results.
%   BLOCKPROCESS(..., DEBUG) Run execution sequentially, i.e. not in parallel, to be
%                       able to step through the code in debug mode.

% Setup
addpath(genpath('common'))
rootfolder = fileparts(mfilename('fullpath'));
poolStartup(poolcores, 'AttachedFiles',{fullfile(rootfolder, 'utils\poolStartup.m')},'debug',debug)
path2data = fullfile(rootfolder,path2data);
writeto   = fullfile(projectpath, 'results');
if ~debug; setupemail; end
fun = func2str(fhandle);
if isstring(varnames), varnames = {varnames}; end

try
    tic
    dd  = dir(fullfile(path2data,'*.mat'));
    if isempty(dd)
        error('No data found in "%s".', path2data)
    end
    N   = numel(dd);
    res = deal(cell(N,1));
    
    % Slice cached
    if isempty(cached)
        cached = cell(N,1);
    elseif ~iscell(cached)
        vnames = setdiff(getVariableNames(cached),'File','stable');
        cached = accumarray(cached.File,(1:size(cached))',[],@(x) {cached(x,vnames)});
    end
    
    % LOOP in parallel
    parfor f = 1:N
        disp(f)
        % Load data
        s     = load(fullfile(path2data,dd(f).name), varnames{:});
        cache = [cached(f,:), f];
        % Convert to tables
        if isfield(s,'data') && isa(s.data,'dataset'), s.data = dataset2table(s.data); end
        if isfield(s,'mst')  && isa(s.mst ,'dataset'), s.mst  = dataset2table(s.mst ); end
        % Apply function
        res{f} = fhandle(s, cache, varargin{:});
    end
    % Collect all results and convert to dataset
    res = cat(1,res{:});
    
    % Export results and notify
    filename = sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),fun);
    save(fullfile(writeto,filename), 'res')
    message  = sprintf('Task ''%s'' terminated in %s',fun,sec2time(toc));
    disp(message)
    if ~debug, sendmail('o.komarov11@imperial.ac.uk', message,''); end
catch err
    filename = fullfile(writeto, sprintf('%s_err.mat',datestr(now,'yyyymmdd_HHMM')));
    save(filename,'err')
    if ~debug
        sendmail('o.komarov11@imperial.ac.uk',sprintf('ERROR in task ''%s''',fun), err.message, {filename})
    end
    rethrow(err)
end
if ~debug
    delete(gcp('nocreate'))
%     rmpref('Internet','SMTP_Password')
end
end