function records = getMaster(symbol, path2mst)
% GETMASTER Loads master records for selected symbols
%
%   GETMASTER(SYMBOL,[PATH2MST])
%       SYMBOL should be a string or a cellstring of symbols
%
%   RECORDS = ...
%       Cell array with a table of the master records for each symbol
%
%       SYMBOL | NAME | CUSIP | ETN | ETA | ETB | ETP | ETX | ETT | ETO| ...
%       cstr   | cstr | cstr  | tf  | tf  | tf  | tf  | tf  | tf  | tf | ...
%
%       ... ETW | ITS | ICODE | SHROUT | UOT | DENOM | TYPE | FDATE
%       ... tf  | tf  | cstr  | double | u32 | char  |  i8  | u32
%
% See also: importMaster

PREFIX_MASTER = 's_';

if isrowchar(symbol)
    symbol = {symbol};
elseif ~iscellstr(symbol)
    error('taq:getMaster:invalid','SYMBOL should be a string or a cellstring.')
end

symbol  = upper(symbol);
files   = fullfile(path2mst,strcat(PREFIX_MASTER, symbol));
n       = numel(symbol);
records = cell(n,1);
for ii = 1:n
    s           = load(files{ii},'-mat');
    records{ii} = s.(char(fieldnames(s)));
end
end
