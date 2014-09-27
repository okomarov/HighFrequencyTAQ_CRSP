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
[~,pos] = ismember(mst(:,{'Id','Date'}), res(:,{'Id','Date'}));
mst.UnID = res.UnID(pos);

% Median and other dailystats
testname = 'dailystats';
try
    res = loadresults(testname);
catch
    res = Analyze(testname,[],mst(:, {'File','Id','Date'}));
end
[~,pos] = ismember(mst(:,{'Id','Date'}), res(:,{'Id','Date'}));
mst.MedPrice = res.MedPrice(pos);
mst.Nrets    = res.Nrets(pos);

% Bad prices days
testname = 'badprices';
try
    res = loadresults(testname);
catch
    res = Analyze(testname,[],mst(:, {'File','Id','Date','MedPrice'}));
end

[~,pos] = ismember(mst(:,{'Id','Date'}), res(:,{'Id','Date'}));
mst.Isbadday = res.Isbadday(pos);

% Bad series
totbad         = accumarray(mst.UnID, mst.Isbadday);
totobs         = accumarray(mst.UnID, mst.To - mst.From +1);
% hist(totbad./totobs,100)
badseries      = totbad./totobs > .1;
badseries(end) = true; % for the unmatched
mst.Isbadday   = mst.Isbadday | badseries(mst.UnID);

% Average time step
testname = 'avgtimestep';
try
    res = loadresults(testname);
catch
    res = Analyze(testname,[], mst(:, {'File','Id','Date','MedPrice'}));
end
[~,pos] = ismember(mst(:,{'Id','Date'}), res(:,{'Id','Date'}));
mst.Timestep = res.Timestep(pos);

% Select on basis of minimum number of observations
% - Worst case 13 trades with an AVERAGE of 30min timestep
ifewtrades   = isnan(res.Timestep) | res.Timestep > 1/48 | mst.Nrets < 12;
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