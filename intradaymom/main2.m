%% Options
OPT_LAGDAY             = 1;
OPT_NOMICRO            = true;
OPT_OUTLIERS_THRESHOLD = 5;
OPT_HASWEIGHTS         = false;

EDGES = serial2hhmmss((9.5:0.5:16)/24);
%% Intraday-average
taq      = loadresults('price_fl');
crsp     = loadresults('dsfquery');
crsp.Prc = abs(crsp.Prc);
cap      = loadresults('cap');

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
s   = loadresults('master');
mst = s.mst;

if OPT_NOMICRO
    idx = isMicrocap(mst,1);
    mst = mst(~idx,:);
end

% Taq open price
taq            = loadresults('price_fl');
[~,pos]        = ismembIdDate(mst.Permno, mst.Date,taq.Permno,taq.Date);
mst.FirstPrice = taq.FirstPrice(pos);

if OPT_HASWEIGHTS
    cap            = getMktCap(mst.Permno,mst.Date,[],[],[],1);
    [idx,pos]      = ismembIdDate(mst.Permno, mst.Date, cap.Permno,cap.Date);
    mst.Cap(idx,1) = cap.Cap(pos(idx));
    mst.Cap(~idx)  = NaN;
end

% Permnos
permnos = unique(mst.Permno);
nseries = numel(permnos);

% Cached
[mst, dates] = cache2cell(mst,  mst.Date);

%%
N   = numel(dates);
avg = NaN(N, numel(EDGES)-1);

tic
poolStartup(8,'AttachedFiles',{'poolStartup.m'})
parfor ii = 2:N
    disp(ii)
    
    % Get 5 min data
    tmp = getTaqData([],[],[],[],[],datapath,mst{ii},false);
    
    % Unstack
    Permno = tmp.Permno(1:79:end);
    HHMMSS = serial2hhmmss(tmp.Datetime(1:79));
    price  = reshape(tmp.Price,79,[]);
    
    % Add first price
    row        = max(sum(isnan(price)),1);
    [~, col]   = ismember(mst{ii}.Permno, Permno);
    pos        = sub2ind(size(price), row(:), col(:));
    price(pos) = mst{ii}.FirstPrice;
    
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
        [~,pos]   = ismember(Permno, mst{ii}.Permno);
        weight    = mst{ii}.Cap(pos);
        weight    = weight(:)'/sum(weight);
        ret       = bsxfun(@times, ret, weight);
        avg(ii,:) = sum(ret,2);
    else
        avg(ii,:) = mean(ret,2);
    end
end
toc

save .\results\avg_tx_30min_vw avg
save .\results\avg_tx_30min_ew avg
%% Plot
avg_vw     = loadresults('avg_tx_30min_vw');
avg_ew     = loadresults('avg_tx_30min_ew');
avg_vw_all = loadresults('avg_ts_vw');
avg_ew_all = loadresults('avg_ts_ew');

% Averages
figure
f = 252*100;

ha = subplot(211);
bar(nanmean(avg_ew)*f)
hold on 
bar(14, mean(avg_ew_all)*f,'r')
title('Average annualized % returns - EW')
set(gca,'XtickLabel',EDGES(1:end-1)/100)

subplot(212)
bar(nanmean(avg_vw)*f)
hold on
bar(14, mean(avg_vw_all)*f,'r')
title('Average annualized % returns - VW')
set(gca,'XtickLabel',EDGES(1:end-1)/100)
legend('half-hour','open-to-close','Location','NorthWest')

% Cumulated returns
figure
dts = yyyymmdd2datetime(dates);
subplot(221)
hl  = plot(dts,cumprod(nan2zero(avg_ew)+1));
title('Cumulative returns - EW')

subplot(222)
sel = [1,2,size(avg_ew,2)];
hl2 = plot(dts,cumprod(nan2zero(avg_ew(:,sel))+1));
set(hl2,{'Color'}, get(hl(sel),'Color'))

subplot(223)
hl = plot(dts, cumprod(nan2zero(avg_vw)+1));
title('Cumulative returns - VW')

subplot(224)
hl2 = plot(dts, cumprod(nan2zero(avg_vw(:,sel))+1));
set(hl2,{'Color'}, get(hl(sel),'Color'))

legend(num2str(EDGES(sel)),'Location','East')