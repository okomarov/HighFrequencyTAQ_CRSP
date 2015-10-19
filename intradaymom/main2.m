%% Options
OPT_LAGDAY             = 1;
OPT_NOMICRO            = true;
OPT_OUTLIERS_THRESHOLD = 5;
OPT_HASWEIGHTS         = false;

EDGES = serial2hhmmss((9.5:0.5:16)/24);
%% Intraday-average
taq  = loadresults('price_fl');
crsp = loadresults('dsfquery');
crsp.Prc = abs(crsp.Prc);
cap  = loadresults('cap');

if OPT_NOMICRO
    idx  = isMicrocap(crsp,1);
    taq  = taq(~idx,:);
    crsp = crsp(~idx,:);
end

% Intersect
[~,ia,ib]   = intersectIdDate(crsp.Permno,crsp.Date, taq.Permno, taq.Date);
crsp        = crsp(ia,:);
taq         = taq(ib,:);
idx         = ismember(xstr2num(cap.Permnos), taq.Permno);
cap.Data    = cap.Data(:,idx);
cap.Permnos = cap.Permnos(idx);
% isequal(crsp.Date, taq.Date)
% isequal(crsp.Permno, taq.Permno)

% Unstack returns
taq.Ret  = taq.LastPrice./taq.FirstPrice-1;
crsp.Ret = crsp.Prc./crsp.Openprc-1;
ret_taq  = sortrows(unstack(taq (:,{'Permno','Date','Ret'}), 'Ret','Permno'),'Date');
ret_crsp = sortrows(unstack(crsp(:,{'Permno','Date','Ret'}), 'Ret','Permno'),'Date');
ret_taq  = ret_taq{:,2:end};
ret_crsp = ret_crsp{:,2:end};

% Filter out outliers
iout           = ret_taq > OPT_OUTLIERS_THRESHOLD;
ret_taq(iout)  = NaN;
ret_crsp(iout) = NaN;

if OPT_HASWEIGHTS
    w        = bsxfun(@rdivide, cap.Data, nansum(cap.Data,2));
    ret_crsp = ret_crsp.*w(1:end-OPT_LAGDAY,:);
    ret_taq  = ret_taq .*w(1:end-OPT_LAGDAY,:);
    avg      = [nansum(ret_crsp,2), nansum(ret_taq,2)];
else
    avg = [nanmean(ret_crsp,2), nanmean(ret_taq,2)];
end
disp(nanmean(avg)*252*100)
%% Data
datapath = '..\data\TAQ\sampled\5min\nobad';

% Index data
master = loadresults('master');

% Permnos
permnos = unique(master.mst.Permno);
nseries = numel(permnos);

% Capitalizations
cap = loadresults('cap');

% NYSE breakpoints
if OPT_NOMICRO
    idx        = isMicrocap(master.mst,1);
    master.mst = master.mst(~idx,:);
end
%% Lag 1 period
w = [NaN(1,nseries); cap.Data(1+OPT_LAGDAY:end,:)];
%% Cache by dates

% master
master.mst     = sortrows(master.mst,'Date','ascend');
[dates,~,subs] = unique(master.mst.Date);
N              = numel(dates);
nrows          = accumarray(subs,1);
mst            = mat2cell(master.mst,nrows,6);

% cap
w = num2cell(w,2);
%%

EDGES = serial2hhmmss((9.5:0.5:16)/24);
avg   = NaN(N, numel(EDGES)-1);

poolStartup(8,'AttachedFiles',{'poolStartup.m'})
tic
parfor ii = 2:N
    disp(ii)
    
    % Get 5 min data
    tmp = getTaqData([],[],[],[],[],datapath,mst{ii},false);
    
    % Unstack
    Permno = tmp.Permno(1:79:end);
    HHMMSS = serial2hhmmss(tmp.Datetime(1:79));
    price  = reshape(tmp.Price,79,[]);
    
    % Filter outliers
    ret      = price(2:end,:)./price(1:end-1,:)-1;
    idx      = abs(ret) > OPT_OUTLIERS_THRESHOLD;
    ret(idx) = NaN;
    
    % Take 30 min returns
    [~,~,row] = histcounts(HHMMSS(1:end-1),EDGES);
    col       = 1:size(ret,2);
    [row,col] = ndgrid(row, col);
    ret       = accumarray([row(:),col(:)], nan2zero(ret(:))+1,[],@prod)-1;

    if OPT_HASWEIGHTS
        % Intersect permnos
        [~,pos]   = ismember(Permno, permnos);
        weight    = w{ii}(pos);
        weight    = weight/sum(weight);
        ret       = bsxfun(@times, ret, weight);
        avg(ii,:) = sum(ret,2);
    else
        avg(ii,:) = mean(ret,2);
    end
end
toc