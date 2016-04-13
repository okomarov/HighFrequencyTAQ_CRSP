function mst = prepareMst(path2data)
% PREPAREMST Loads master, adds Permno and median price

% Load big master file
commonroot = fullfile(fileparts(mfilename('fullpath')),'..'); 
commonres  = fullfile(commonroot, 'results');
if nargin < 1 || isempty(path2data)
    path2data  = fullfile(commonroot, 'data\TAQ');
end
load(fullfile(path2data,'master'),'-mat')

mst = addPermno(mst);

% Median price
testname = 'medianprice';
try
    res = loadresults(testname, commonres);
catch
    res = Analyze(testname,[],mst(:, {'File','Id','Date'}));
end
[~,pos]      = ismembIdDate(mst.Id, mst.Date, res.Id, res.Date);
mst.MedPrice = res.MedPrice(pos);

% Drop incomplete days
mst.Isbadday = isprobdate(mst.Date);
end