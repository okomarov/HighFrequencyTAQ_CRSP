%% Check Open/Close in TAQ vs CRSP
OPT_NOMICRO = true;
OPT_OUTLIERS_THRESHOLD = 1;
OPT_LAGDAY = 1;
taq = loadresults('sampleFirstLast','..\results');
% taq = load('D:\TAQ\HF\results\20150921_0044_sampleFirstLast.mat');
% taq = taq.res;

if OPT_NOMICRO
    idx = isMicrocap(taq,'LastPrice',OPT_LAGDAY);
    taq = taq(~idx,:);
end

crsp      = loadresults('dsfquery','..\results');
crsp      = crsp(crsp.Prc > 0,:);
[~,ia,ib] = intersectIdDate(crsp.Permno,crsp.Date, taq.Permno, taq.Date);
crsp      = crsp(ia,:);
taq       = taq(ib,:);

% Filter outliers
taq.TAQret   = taq.LastPrice./taq.FirstPrice-1;
iout         = taq.TAQret          > OPT_OUTLIERS_THRESHOLD |...
               1./(taq.TAQret+1)-1 > OPT_OUTLIERS_THRESHOLD;
taq(iout,:)  = [];
crsp(iout,:) = [];

% Comparison table
cmp         = [crsp(:,{'Date','Permno','Openprc'}) taq(:,'FirstPrice'), ...
               crsp(:,{'Bid','Ask','Prc'}),taq(:,{'LastPrice','TAQret'})];
cmp.CRSPret = cmp.Prc./cmp.Openprc-1;

% Filter NaNs
inan = isnan(cmp.TAQret) | isnan(cmp.CRSPret);
cmp  = cmp(~inan,:);

corr(cmp.TAQret, cmp.CRSPret)
ir = mean(rets(:,1)-rets(:,2))/std(rets(:,1)-rets(:,2));

[dts,~,subs] = unique(cmp.Date);
rets = [accumarray(subs,cmp.TAQret,[],@mean), accumarray(subs,cmp.CRSPret,[],@mean)];
plot(yyyymmdd2datetime(dts), cumprod(rets+1))

retdiff = abs(cmp.CRSPret - cmp.TAQret);
idx     = retdiff > eps*1e12;
boxplot(retdiff(idx))

%% Check new/old ff49
load('D:\TAQ\HF\intradaymom\results\bck\FF49.mat')

taq  = loadresults('price_fl');
ff49 = getFF49IndustryCodes(taq,1);
ff49 = struct('Permno', xstr2num(getVariableNames(ff49(:,2:end))), ...
    'Dates', ff49{:,1},...
    'Data', ff49{:,2:end});

% Intersect permno
[~,ia,ib] = intersect(industry.Permno, ff49.Permno);
industry.Data = industry.Data(:,ia);
industry.Permno = industry.Permno(ia);
ff49.Data = ff49.Data(:,ib);
ff49.Permno = ff49.Permno(ib);

% intersect date
[~,ia,ib] = intersect(industry.Date, ff49.Dates);
industry.Data = industry.Data(ia,:);
industry.Date = industry.Date(ia);
ff49.Data = ff49.Data(ib,:);
ff49.Dates = ff49.Dates(ib);

[r,c] = find(ff49.Data ~= industry.Data & ff49.Data ~= 0)
%% CRSP ind correlations
OPT_LAGDAY = 1;
OPT_HASWEIGHTS = true;

crsp     = loadresults('dsfquery');
crsp.Prc = abs(crsp.Prc);
idx      = isMicrocap(crsp,'Prc',OPT_LAGDAY);
crsp     = crsp(~idx,:);

% Get market caps
cap = getMktCap(crsp,OPT_LAGDAY,true);
cap = struct('Permnos', {getVariableNames(cap(:,2:end))}, ...
    'Dates', cap{:,1},...
    'Data', cap{:,2:end});

% FF49-industries classification
ff49 = getFF49IndustryCodes(crsp,1);
ff49 = struct('Permnos', {getVariableNames(ff49(:,2:end))}, ...
    'Dates', ff49{:,1},...
    'Data', ff49{:,2:end});

% Unstack returns
crsp.Ret = crsp.Prc./crsp.Openprc-1;
ret_crsp = sortrows(unstack(crsp(:,{'Permno','Date','Ret'}), 'Ret','Permno'),'Date');
ret_crsp = ret_crsp{:,2:end};

if OPT_HASWEIGHTS
    w = bsxfun(@rdivide, cap.Data, nansum(cap.Data,2));
else
    w = repmat(1./sum(~isnan(ret_taq),2), 1,size(ret_taq,2));
end
ret_crsp_w = ret_crsp.*w;

nobs    = numel(ff49.Dates);
nseries = numel(ff49.Permnos);
row     = repmat((1:nobs)', 1, nseries);
subs    = [row(:), double(ff49.Data(:)+1)];

ptfret_vw = accumarray(subs, ret_crsp_w(:),[],@nansum);
c = corr(ptfret_vw); c(logical(eye(50))) = NaN