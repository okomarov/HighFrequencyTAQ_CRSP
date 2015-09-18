%% Options
OPT_HASWEIGHTS       = false;
OPT_LAGDAY           = 1;
OPT_PTFNUM           = 10;
OPT_NOMICRO          = false;
OPT_EDGES_BAD_PRICES = [.5,1.5];

%% Data
datapath = '..\data\TAQ\sampled\5min\new';

% Index data
master = loadresults('master');

% First and last price
price_fl = loadresults('price_fl');

% Permnos
permnos = unique(master.mst.Permno);
nseries = numel(permnos);

% Capitalizations
if OPT_HASWEIGHTS
    % First and last price
    cap = loadresults('cap');
end

% NYSE breakpoints
if OPT_NOMICRO
    bpoints = loadresults('ME_breakpoints_TXT','..\results');
    idx     = ismember(bpoints.Date, unique(master.mst.Date/100));
    bpoints = bpoints(idx,{'Date','Var3'});
end
%% Lag 1 period
if OPT_HASWEIGHTS
    w = [NaN(1,nseries); cap.Data(1+OPT_LAGDAY:end,:)];
end
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
price_fl = mat2cell(price_fl,nrows,6);

% cap
if OPT_HASWEIGHTS
    w = num2cell(w,2);
else
    w = cell(N,1);
end
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
    
    % Signal: Filled back half-day ret
    past_ret = flipud(nanfillts(price(37:-1:1,:)));
    past_ret = past_ret(end,:)./past_ret(1,:)-1;
    
    % Filter out microcaps
    
    % hpr with 5 min skip
    ret       = price(end,:)./price(38,:)-1;
    ptf(ii,:) = portfolio_sort(past_ret,ret, 'PortfolioNumber',10,'Weights',w{ii});
end
toc
stratstats(ptf)