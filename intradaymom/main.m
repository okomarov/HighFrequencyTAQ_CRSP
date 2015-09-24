%% Options
OPT_HASWEIGHTS    = true;
OPT_LAGDAY        = 1;
OPT_PTFNUM        = 10;
OPT_PTFNUM_DOUBLE = [5,5];
OPT_INDEP_SORT    = false;
OPT_NOMICRO       = true;

%% Data
datapath = '..\data\TAQ\sampled\5min\nobad';

% Index data
master = loadresults('master');

% First and last price
price_fl = loadresults('price_fl');

% Permnos
permnos = unique(master.mst.Permno);
nseries = numel(permnos);

% Capitalizations
cap = loadresults('cap');

% NYSE breakpoints
if OPT_NOMICRO
    bpoints = loadresults('ME_breakpoints_TXT','..\results');
    idx     = ismember(bpoints.Date, unique(master.mst.Date/100));
    bpoints = bpoints(idx,{'Date','Var3'});
end
%% Lag 1 period
w = [NaN(1,nseries); cap.Data(1+OPT_LAGDAY:end,:)];

if OPT_NOMICRO
    bpoints.Var3 = [NaN(OPT_LAGDAY,1); bpoints.Var3(1:end-OPT_LAGDAY)];
end
%% Cache by dates

% master
master.mst     = sortrows(master.mst,'Date','ascend');
[dates,~,subs] = unique(master.mst.Date);
N              = numel(dates);
nrows          = accumarray(subs,1);
mst            = mat2cell(master.mst,nrows,6);

% price first last
price_fl       = sortrows(price_fl,'Date','ascend');
[dates,~,subs] = unique(price_fl.Date);
nrows          = accumarray(subs,1);
price_fl       = mat2cell(price_fl,nrows,size(price_fl,2));

% cap
w = num2cell(w,2);
%%

ptf = NaN(N,OPT_PTFNUM);

poolStartup(8,'AttachedFiles',{'poolStartup.m'})
tic
parfor ii = 2:N
    disp(ii)
    
    % Get sampled data
    tmp           = getTaqData([],[],[],[],[],datapath,mst{ii},false);
    % unstack
    price         = NaN(79, nseries);
    [~,col]       = ismember(tmp.Permno, permnos);
    [un,irow,row] = unique(serial2hhmmss(tmp.Datetime));
    pos           = (col-1)*79+row;
    price(pos)    = tmp.Price;
        
    price_first      = NaN(1,nseries);
    [idx,col]        = ismember(price_fl{ii}.Permno, permnos);
    price_first(col) = price_fl{ii}.FirstPrice(idx);
    
    price_last      = NaN(1,nseries);
    price_last(col) = price_fl{ii}.LastPrice(idx);
    
    % Signal: Filled back half-day ret
    past_ret = price(FORMATION,:)./price_first-1;

    % hpr with 5 min skip
    hpr = price_last./price(FORMATION+2,:)-1;
    
    % Filter microcaps
    if OPT_NOMICRO
        nyseCap  = bpoints.Var3(ismember(bpoints.Date, price_fl{ii}.Date/100));
        idx      = price_first < 5 | w{ii} < nyseCap;
        hpr(idx) = NaN;
    end
    
    if OPT_HASWEIGHTS
        weight = w{ii};
    else
        weight = [];
    end
    % PTF ret
    ptf(ii,:) = portfolio_sort(hpr,past_ret, 'PortfolioNumber',OPT_PTFNUM, 'Weights',weight);
    
    % PTF ret
    [ptf2(ii,:), bin2(ii,:)] = portfolio_sort(hpr,{w{ii},past_ret}, 'PortfolioNumber',OPT_PTFNUM_DOUBLE,...
        'Weights',weight,'IndependentSort',OPT_INDEP_SORT);
end
toc

t = stratstats(dates, [ptf, ptf(:,1)-ptf(:,end)] ,'d',0);
t{:,:}'

t2 = stratstats(dates, ptf2 ,'d',0);
t2{:,:}'
reshape(t2.Annret, OPT_PTFNUM_DOUBLE)'
% t.Properties.VariableNames'

OPT_HASWEIGHTS
