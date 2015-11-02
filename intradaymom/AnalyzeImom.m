function [res, filename] = AnalyzeImom(fun, varnames, cached, path2data, debug, poolcores, varargin)
% ANALYZE Executes specified fun in parallel on the whole database (all .mat files)
%
%   ANALYZE(FUN, VARNAMES) FUN should a string with the name of one of
%                          the following sub-functions:
%                               - 'dailystats'
%                               - 'badprices'
%                               - 'avgtimestep'
%                          VARNAMES is a cell-array of strings (or string)
%                          with the VarNames of the dataset with the results
%
%   ANALYZE(..., PATH2DATA) If you wanna use other than '.\data\TAQ\'
%                           files (default), then specify a different
%                           PATH2DATA, e.g '.\data\TAQ\sampled\5min'
%
%   ANALYZE(..., CACHED) Some FUN might require pre-cached results which where
%                        run on the whole database.
%                        Check the specific sub-function for the format of the
%                        needed CACHED results.
%   ANALYZE(..., DEBUG) Run execution sequentially, i.e. not in parallel, to be
%                       able to step through the code in debug mode.
if nargin < 2 || isempty(varnames),  varnames  = '';             end
if nargin < 3,                       cached    = [];             end
if nargin < 4 || isempty(path2data); path2data = '.\data\TAQ\';  end
if nargin < 5 || isempty(debug);     debug     = false;          end
if nargin < 6 || isempty(poolcores); poolcores = 8;              end

fhandles = {@getPrices};
    
[hasFunc, pos] = ismember(fun, cellfun(@func2str,fhandles,'un',0));
if ~hasFunc
    error('Unrecognized function "%s".', fun)
end
fun             = fhandles{pos};
projectpath     = fileparts(mfilename('fullpath'));
[res, filename] = blockprocess(fun,projectpath, varnames, cached,path2data,debug,poolcores,varargin{:});
end

%% Subfunctions
