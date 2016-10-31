function record = getTaqMaster(idtype,id,date,ids)

if nargin < 4
    try
        master = load(fullfile('data\TAQ','master'),'-mat','ids');
    catch
        master = load(fullfile('..\data\TAQ','master'),'-mat','ids');
    end
    ids = master.ids;
end
try
    mstp = loadresults('masterPermno');
catch
    mstp = loadresults('masterPermno','..\results');
end
if strcmpi(idtype,'permno')
    permno = id;
    id     = unique(mstp.Id(mstp.Permno == id));
elseif strcmpi(idtype,'id')
    permno = mstp.Permno(mstp.Id == id & mstp.Date == date);
end

try
    master = loadresults('TAQmaster');
catch
    master = loadresults('TAQmaster','..\results');
end

% TAQ
tmp = master(ismember(master.SYMBOL,ids{id}),:);
tmp = sortrows(tmp, 'FDATE');
if size(tmp,1) == 1
    record = tmp(:,1:3);
else
    [~,~,bin] = histcounts(date,[double(tmp.FDATE)-0.5; inf]);
    record       = tmp(bin,1:3);
end
record.Date = date;
record.Id   = id;
record.Permno = permno;
end
