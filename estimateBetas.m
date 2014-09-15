function betas = estimateBetas(lookback, freq, useproxy)
if nargin < 1 || isempty(lookback), lookback = 1;     end
if nargin < 2 || isempty(freq),     freq     = 5;     end
if nargin < 3 || isempty(useproxy), useproxy = false; end

% Sample if data doesn't exist
path2data = sprintf('.\\data\\TAQ\\sampled\\%dmin', freq);
fmtname   = sprintf('S%dm_%%04d.mat',freq);
if exist(path2data,'dir') ~= 7 || isempty(ls(path2data))
    step    = freq/(60*24);
    grid    = (9.5/24:step:16/24)';
    sampleData(grid, path2data, fmtname);
end 

% Use self-built sp500 proxy
if useproxy
    name = sprintf('sp500proxy%dm',freq);
    try
        loadresults(name, 'sp500')
    catch
        sp500 = sp500intraday(path2data);
    end
% Use spyders    
else
    name = sprintf('spysampled%dm',freq);
    try
        loadresults(name, 'sp500')
    catch
        master = load(fullfile(path2data,'master'),'-mat');
        sp500  = getTaqData(master, 'SPY',[],[],'Price',path2data);
        save(sprintf('.\results\%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),name), 'sp500')
    end
end

    % SP500 ret (zeroing overnight)
    spret = [sp500.Datetime(2:end) sp500.Price(2:end)./sp500.Price(1:end-1)-1];
    spret = spret(diff(rem(sp500.Datetime,1)) >= 0,:);
    
    % Cache SP500 returns by days
    load(fullfile(path2data,'master'),'-mat','mst');
    nfiles = max(mst.File);
    cached = cell(nfiles,1);
    
    dates           = fix(spret(:,1));
    [spdays,~,subs] = unique(dates,'stable');
    spret           = mat2cell(spret(:,2), accumarray(subs,1),1);
    
    unMstDates = accumarray(mst.File, mst.Date,[],@(x){yyyymmdd2serial(unique(x))});
    for ii = 1:nfiles
        pos        = ismembc2(unMstDates{ii}, spdays);
        nnzero     = pos ~= 0;
        isp        = ismembc(spdays, unMstDates{ii});
        cached{ii} = {spret(isp) spdays(pos(nnzero))};
    end
    
    % Calculate betas
    bcomp = Analyze('betacomponents', [], cached, fullfile(path2data,fmtname));
end