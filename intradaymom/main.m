%% Options
OPT_HASWEIGHTS       = false;
OPT_LAGDAY           = 1;
OPT_PTFNUM           = 10;
OPT_NOMICRO          = false;
OPT_EDGES_BAD_PRICES = [.5,1.5]; 

%% Data
datapath = '..\data\TAQ\sampled\5min\new';

% Sample first and last price
mst = selectAndFilterTrades(OPT_EDGES_BAD_PRICES);
if isempty(OPT_EDGES_BAD_PRICES)
    mst = mst(:, {'File','Id','Permno','Date','Isbadday','Isfewobs'});
else
    mst = mst(:, {'File','Id','Permno','Date','MedPrice', 'Isbadday','Isfewobs'});
end
Analyze('sampleFirstLast',[],mst,[],1, struct('edges',OPT_EDGES_BAD_PRICES))

% Index data
master     = load(fullfile(datapath,'master'),'-mat');
shrcd      = loadresults('shrcd','../results');
idx        = iscommonshare(master.mst, shrcd);
master.mst = master.mst(idx,:);

% Permnos
permnos = unique(master.mst.Permno);
nseries = numel(permnos);

% Capitalizations
if OPT_HASWEIGHTS
    cap = getMktCap(master.mst.Permno, master.mst.Date,false,false,1);
    cap = struct('Permnos', {getVariableNames(cap(:,2:end))}, ...
        'Dates', cap{:,1},...
        'Data', cap{:,2:end});
end

% NYSE breakpoints
if OPT_NOMICRO
    try
        bpoints = loadresults('ME_breakpoints','..\results');
    catch
        bpoints = importFrenchData('ME_Breakpoints_TXT.zip','..\results');
    end
end
%% Lag 1 day
if OPT_HASWEIGHTS
    w = [NaN(1,nseries); cap.Data(1+OPT_LAGDAY:end,:)];
end
%% Cache by dates

% master
master.mst     = sortrows(master.mst,'Date','ascend');
[dates,~,subs] = unique(master.mst.Date);
N              = numel(dates);
nrows          = accumarray(subs,1);
mst            = mat2cell(master.mst,nrows,6);

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
    
    % Get data
    tmp           = getTaqData([],[],[],[],[],datapath,mst{ii},false);
    % unstack
    slice         = NaN(79, nseries);
    [~,col]       = ismember(tmp.Permno, permnos);
    [un,irow,row] = unique(serial2hhmmss(tmp.Datetime));
    pos           = (col-1)*79+row;
    slice(pos)    = tmp.Price;
    
    % Signal: Filled back half-day ret
    past_ret = flipud(nanfillts(slice(37:-1:1,:)));
    past_ret = past_ret(end,:)./past_ret(1,:)-1;
    
    % hpr with 5 min skip
    ret       = slice(end,:)./slice(38,:)-1;
    ptf(ii,:) = portfolio_sort(past_ret,ret, 'PortfolioNumber',10,'Weights',w{ii});
end
toc
stratstats(ptf)