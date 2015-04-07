function sampleData(grid, writeto, fmtname)

% Write path handling
[isok, msg, msgid] = mkdir(writeto);
direxists = strcmpi('MATLAB:MKDIR:DirectoryExists',msgid);
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
path2data = '.\data\TAQ';
load(fullfile(path2data,'master'),'-mat')

% Map unique ID to mst
testname = 'uniqueID';
try
    res = loadresults(testname);
catch
    res = mapUnid2mst(mst, ids);
end
[~,pos] = ismembIdDate(mst.Id, mst.Date, res.Id, res.Date);
mst.UnID = res.UnID(pos);

% Median price
testname = 'medianprice';
try
    res = loadresults(testname);
catch
    res = Analyze(testname,[],mst(:, {'File','Id','Date'}));
end
[~,pos]      = ismembIdDate(mst.Id, mst.Date, res.Id, res.Date);
mst.MedPrice = res.MedPrice(pos);

% Bad prices days
testname = 'badprices';
try
    res = loadresults(testname);
catch
    dailycut = 0.5;
    res = Analyze(testname,[],mst(:, {'File','Id','Date','MedPrice'}),[],[],dailycut);
end
[~,pos]      = ismembIdDate(mst.Id, mst.Date, res.Id, res.Date);
mst.Isbadday = res.Isbadday(pos);
mst.Nbadtot  = res.Nbadtot(pos);

% Bad series
nobs           = mst.To - mst.From +1;
totbad         = accumarray(mst.UnID(mst.Isbadday), nobs(mst.Isbadday),[max(mst.UnID),1]);
totobs         = accumarray(mst.UnID, mst.To - mst.From +1);
badseries      = totbad./totobs > .1;
badseries(end) = true; % for the unmatched
mst.Isbadday   = mst.Isbadday | badseries(mst.UnID);

% Count losing obs with timestamp consolidation
testname = 'consolidationcounts';
try
    res = loadresults(testname);
catch
    res = Analyze(testname,[],mst(:, {'File','Id','Date','MedPrice','Isbadday'}));
end
[~,pos]           = ismembIdDate(mst.Id, mst.Date, res.Id, res.Date);
mst.Nconsolidated = res.Nconsolidated(pos);

% Count number of time buckets in a day that have a trade
testname = 'NumTimeBuckets';
try
    res = loadresults(testname);
catch
    res = Analyze('NumTimeBuckets',[],mst(:, {'File','Id','Date','MedPrice','Isbadday'}));
end
[~,pos]            = ismembIdDate(mst.Id, mst.Date, res.Id, res.Date);
mst.NumTimeBuckets = res.NumTimeBuckets(pos);

% Select with minimum number of observations
% ngoodtrades  = nobs - mst.Nbadtot - mst.Nconsolidated;
% ifewtrades   = ngoodtrades < 15;
ifewtrades   = mst.NumTimeBuckets < 7;
perfew       = accumarray(mst.UnID, ifewtrades)./accumarray(mst.UnID, 1) > .5;
mst.Timestep = ifewtrades | perfew(mst.UnID);

%% Sample at x min
fprintf('%s: sampling.\n', mfilename)
mst = mst(:, {'File','Id', 'Date', 'UnID','MedPrice','Isbadday','Timestep'});
opt = struct('grid',grid, 'writeto', writeto, 'fmtname', fmtname);
clearvars -except mst opt
Analyze('sample',[], mst,[],[],opt);

% Sort mst by from
fprintf('%s: sorting all mst by ''From''.\n', mfilename)
sortmst(opt.writeto)

% Make master file
fprintf('%s: making master records file.\n', mfilename)
makeMasterFile(opt.writeto)
end