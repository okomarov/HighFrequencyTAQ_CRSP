function mst = prepareMst()

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

% Median price
testname = 'medianprice';
try
    res = loadresults(testname, commonres);
catch
    res = Analyze(testname,[],mst(:, {'File','Id','Date'}));
end
[~,pos]      = ismembIdDate(mst.Id, mst.Date, res.Id, res.Date);
mst.MedPrice = res.MedPrice(pos);
end