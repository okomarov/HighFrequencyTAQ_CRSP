function sampleData(grid, writeto, fmtname)

% Write path handling
[status, msg, msgid] = mkdir(writeto);
if ~status || strcmpi('MATLAB:MKDIR:DirectoryExists',msgid)
    if ~isempty(ls(writeto))
        msg = 'the path ''%s'' already exists and is not empty.\n\t\t Data might get overwritten.';
        warning('sampleData:nonEmptyWriteDir',msg,writeto)
    else
        error(msgid,msg)
    end
end
%% Selection/filtering

% Load big master file
path2data = '.\data\TAQ';
load(fullfile(path2data,'master'),'-mat')

% Map unique ID to mst
testname = 'uniqueID';
try
    loadresults(testname,'res')
catch
    res = mapUnid2mst(mst, ids);
end
[~,pos] = ismember(mst(:,{'Id','Date'}), res(:,{'Id','Date'}));
mst     = [mst, res(pos, 'UnID')];

% Median and other dailystats
testname = 'dailystats';
try
    loadresults(testname,'res')
catch
    res = Analyze(testname,{'Min','Max','MedPrice','Nrets'});
end
if isa(res, 'dataset'), res = dataset2table(res); end 
[~,pos] = ismember(mst(:,{'Id','Date'}), res(:,{'Id','Date'}));
mst     = [mst, res(pos,{'MedPrice','Nrets'})];

% Bad prices days
testname = 'badprices';
try
    loadresults(testname,'res')
catch
    res = Analyze(testname,'Baddays',mst(:, {'File','MedPrice'}));
end
if isa(res, 'dataset'), res = dataset2table(res); end
[~,pos] = ismember(mst(:,{'Id','Date'}), res(:,{'Id','Date'}));
mst     = [mst, res(pos,'Baddays')];

% Bad series
totbad         = accumarray(mst.UnID, mst.Baddays);
totobs         = accumarray(mst.UnID, mst.To - mst.From +1);
% hist(totbad./totobs,100)
badseries      = totbad./totobs > .1;
badseries(end) = true; % for the unmatched
mst.Baddays    = mst.Baddays | badseries(mst.UnID);

% Average time step
testname = 'avgtimestep';
try
    loadresults(testname,'res')
catch
    res = Analyze(testname,'Timestep', mst(:, {'File','MedPrice'}));
end
if isa(res, 'dataset'), res = dataset2table(res); end
[~,pos] = ismember(mst(:,{'Id','Date'}), res(:,{'Id','Date'}));
mst     = [mst, res(pos,'Timestep')];

% Select on basis of minimum number of observations
% - Worst case 13 trades with an AVERAGE of 30min timestep
ifewtrades   = isnan(res.Timestep) | res.Timestep > 1/48 | mst.Nrets < 12;
perfew       = accumarray(mst.UnID, ifewtrades)./accumarray(mst.UnID, 1) > .5;
mst.Timestep = ifewtrades | perfew(mst.UnID);

%% Sample at x min
mst = mst(:, {'File','Id','UnID','MedPrice','Baddays','Timestep'});
opt = struct('grid',grid, 'writeto', writeto, 'fmtname', fmtname);
Analyze('sample',[], mst,[],1,opt);

% Sort mst by from
sortmst(opt.writeto)

% Make master file
makeMasterFile(opt.writeto)

end