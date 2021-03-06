function sampleData(grid, writeto, fmtname, badPriceMult, consolidateType)
% sampleData(grid, writeto, fmtname, edgesBadPrices)

% Sampling params
if nargin < 1 || isempty(grid),    grid            = (9.5/24:5/(60*24):16/24)'; end
if nargin < 2 || isempty(writeto), writeto         = '.\data\TAQ\sampled\5min'; end
if nargin < 3 || isempty(fmtname), fmtname         = 'S5m_%04d.mat';            end
if nargin < 4,                     badPriceMult    = [];                        end
if nargin < 5,                     consolidateType = 'volumeWeighted';          end

% Write path handling
[isok, msg, msgid] = mkdir(writeto);
direxists          = strcmpi('MATLAB:MKDIR:DirectoryExists',msgid);
if ~isok
    error(msgid,msg)
elseif direxists 
    if numel(dir(writeto)) > 2
        msg = 'the path ''%s'' already exists and is not empty.\n\t\t Data might get overwritten.';
        warning('sampleData:nonEmptyWriteDir',msg,writeto)
    end
else
    fprintf('%s: created new dir ''%s''.\n', mfilename, writeto)
end
%% Selection/filtering
fprintf('%s: selection and filtering.\n', mfilename)

% Load big master file
commonroot = fullfile(fileparts(mfilename('fullpath')),'..'); 
commonres  = fullfile(commonroot, 'results');
path2data  = fullfile(commonroot, 'data\TAQ');
load(fullfile(path2data,'master'),'-mat')

% Map unique ID to mst
testname = 'masterPermno';
try
    res = loadresults(testname,commonres);
catch
    res = mapPermno2master;
end
[~,pos]    = ismembIdDate(mst.Id, mst.Date, res.Id, res.Date);
mst.Permno = res.Permno(pos);

% Nobs
mst.Nobs = mst.To - mst.From +1;

if ~isempty(badPriceMult)
    % Median price
    testname = 'medianprice';
    try
        res = loadresults(testname, commonres);
    catch
        res = Analyze(testname,[],mst(:, {'File','Id','Date'}));
    end
    [~,pos]      = ismembIdDate(mst.Id, mst.Date, res.Id, res.Date);
    mst.MedPrice = res.MedPrice(pos);
    keepflds     = {'File','Id','Date','MedPrice'};
else
    keepflds = {'File','Id','Date'};
end

% Drop incomplete days
mst.Isbadday = isprobdate(mst.Date);

%% Sample at x min
fprintf('%s: sampling.\n', mfilename)
if ~isempty(badPriceMult)
    mst = mst(:, {'File','Id','Permno','Date','MedPrice','Isbadday'});
else
    mst = mst(:, {'File','Id','Permno','Date','Isbadday'});
end
opt = struct('grid',grid, 'writeto', writeto, 'fmtname', fmtname,...
    'BadPriceMultiplier',badPriceMult,'TimestampConsolidation',consolidateType);
clearvars -except mst opt
Analyze('sample',[], mst,[],[],[],opt);

% Sort mst by from
fprintf('%s: sorting all mst by ''From''.\n', mfilename)
sortmst(opt.writeto)

% Make master file
fprintf('%s: making master records file.\n', mfilename)
makeMasterFile(opt.writeto)
end