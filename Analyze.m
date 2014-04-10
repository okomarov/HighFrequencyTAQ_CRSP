function Analyze(fun, varnames, cached, debug)
if nargin < 3; cached = [];    end
if nargin < 4; debug  = false; end   

% Simply call the specific subroutine
cd C:\HFbetas
addpath(genpath('.\utils'))
writeto = '.\results\';

% Open matlabpool
if matlabpool('size') == 0 && ~debug
     matlabpool('open', 4, 'AttachedFiles',{'.\utils\poolStartup.m'})
end

% Get email credentials if not in debug
if ~debug; setupemail; end
    
fhandle = str2func(fun);

try
    tic
    d   = '.\data\TAQ';
    dd  = dir(fullfile(d,'*.mat'));
    N   = numel(dd);
    res = deal(cell(N,1));
    
    % Slice cached
    if isempty(cached)
        cached = cell(N,1);
    elseif ~iscell(cached)
        cached = accumarray(cached.File,(1:size(cached))',[],@(x) {cached(x,2:end)});
    end
    
    % LOOP in parallel
    parfor f = 1:N
        disp(f)
        % Load data
        s      = load(fullfile(d,dd(f).name));
        % Apply function
        res{f} = fhandle(s, cached{f});
    end
    % Collect all results and convert to dataset
    res = cat(1,res{:});
    if ~isempty(varnames)
        res = mat2dataset(res,'VarNames',varnames);
    end
    
    % Export results and notify
    save(fullfile(writeto,sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),fun)), 'res')
    message = sprintf('Task ''%s'' terminated in %s',fun,sec2time(toc));
    disp(message)
    if ~debug, sendmail('o.komarov11@imperial.ac.uk', message,''); end
catch err
    filename = fullfile(writeto, sprintf('%s_err.mat',datestr(now,'yyyymmdd_HHMM')));
    save(filename,'err')
    if ~debug
        sendmail('o.komarov11@imperial.ac.uk',sprintf('ERROR in task ''%s''',fun), err.message, {filename})
    end
    rethrow(err)
end
if ~debug
    matlabpool close
    rmpref('Internet','SMTP_Password')
end
end
%% Count bad prices
function res = badprices(s, cached)

% Detect bad prices
% If number of bad prices of that day is 50% then is a bad day
% If the number of bad days per series is more than 50% then bad series

dailycut = .5;

% Number of observations per day
nobs = double(s.mst.To - s.mst.From + 1);

% STEP 1) Selection
inan = selecttrades(s.data);

% STEP 2) Bad prices are 1.5x or .65x the daily median (net of selection nans)
[~,ibadprice] = histc(s.data.Price./RunLength(cached.MedPrice,nobs), [.65,1.51]);
ibadprice     = ibadprice ~= 1;

% STEP 3) Bad days
 res = accumarray(RunLength((1:size(s.mst,1))',nobs), inan | ibadprice) > ceil(dailycut*nobs);
end
%% Check how many bad prices
function res = cleanprices(s,cached)
% Number of observations per day
nobs = double(s.mst.To - s.mst.From + 1);

% STEP 1) Selection
inan   = selecttrades(s.data);
selper = nnz(inan)/numel(s.data.Price)*100;
% fprintf('Selection: %.1f\n',selper)

% STEP 2) Bad prices prices are 1.5x or .65x the daily median (net of selection nans)
[~,igoodprice] = histc(s.data.Price./RunLength(cached.MedPrice,nobs), [.65,1.51]);
inan           = inan | igoodprice ~= 1;
% fprintf('Bad prices: %.1f\n',nnz(inan)/numel(s.data.Price)*100)

% STEP 3) Bad series/days
inan      = inan | RunLength(cached.Baddays,nobs);
badserper = nnz(inan)/numel(s.data.Price)*100;
% fprintf('Bad series/days: %.1f\n',badserper)

% STEP 4) Filter out days with < 30min avg timestep or securities with 50% fewtrades days
inan = inan | RunLength(cached.Timestep,nobs);
% fprintf('Timestep: %.1f\n',nnz(inan)/numel(s.data.Price)*100)

% STEP 5) Clean prices
s.data.Price(inan) = NaN;

res = NaN(1,2,'single');
if ~all(inan)
    % Returns
    ret  = s.data.Price(2:end)./s.data.Price(1:end-1)-1;
    % Excluding overnight
    ret(setdiff(s.mst.From-1,[-1,0])) = NaN;
    % MinMax
    [res(1),res(2)] = MinMaxElem(ret);
end
fprintf('\t\t\t\t\t\t\t %.1f - %.1f\n',selper,badserper)
end
%% Betas
function res = betas(s,cached)
% Number of observations per day
nobs = double(s.mst.To - s.mst.From + 1);
nmst = size(s.mst,1);
sp500 = cached{2};
cached = cached{1};
% STEP 1) Selection
inan = selecttrades(s.data);

% STEP 2) Bad prices prices are 1.5x or .65x the daily median (net of selection nans)
[~,iprice] = histc(s.data.Price./RunLength(cached.MedPrice,nobs), [.65,1.51]);
iprice     = iprice ~= 1;
inan       = inan | iprice;

% STEP 3) Bad series/days
inan = inan | RunLength(cached.Baddays,nobs);

% STEP 4) Filter out days with < 30min avg timestep or securities with 50% fewtrades days
inan = inan | RunLength(cached.Timestep,nobs);

% STEP 5) Clean prices
s.data.Price(inan) = NaN;

res = NaN(nmst,1,'single');
if ~all(inan)

    % STEP 6) Median prices for same timestamps
    mstrow           = RunLength((1:nmst)',nobs);
    mstrow           = mstrow(~inan);
    [unTimes,~,subs] = unique(mstrow + hhmmssmat2serial(s.data.Time(~inan,:)));
    price            = accumarray(subs, s.data.Price(~inan),[],@fast_median);
        
    % STEP 7) Sample on fixed grid (easier to match sp500) 
    grid  = (9.5/24:5/(60*24):16/24)';
    ngrid = numel(grid);
    
%     biggrid = bsxfun(@plus, grid, unique(fix(unTimes))');
%     [price,pricedt] = realized_price_filter(price,unTimes,'Unit','Fixed',biggrid(:))
  
    % Sample price and 
    price = fixedsampling(double([unTimes, price]), 'Previous', grid);
    
    % Collect sp500 into cell per day array
    sp500date2cell = fix(sp500.Datetime(ngrid:ngrid:end));
    sp500          = mat2cell(sp500.Price,repmat(ngrid,size(sp500,1)/ngrid,1),1);
    [~,pos]        = ismember(yyyymmdd2serial(double(s.mst.Date(fix(price(ngrid:ngrid:end,1))))), sp500date2cell);
    sp500          = cat(1,sp500{pos});
        
    % STEP 8) Sample SP500
%     sp500 = fltprice(double(sp500),'FixedTime', 'Previous', grid);
%     biggrid2 = bsxfun(@plus, grid, unique(fix(sp500.Datetime))');
%     [sp500,sp500dt] = realized_price_filter(sp500.Price,sp500.Datetime,'Unit','Fixed',biggrid2(:))
   
    % Returns
    ret      = price(2:end,2)./price(1:end-1,2)-1;
    retSP500 = sp500(2:end)./sp500(1:end-1)-1;
    nret     = numel(ret);
    prodret  = ret .* retSP500;
    % Get rid of overnight and NaNs
    subs     = reshape(repmat(1:numel(pos),ngrid,1),[],1);
    ikeep    = ~isnan(prodret);
    ikeep(ngrid:ngrid:nret) = false;
    % Realised beta
    res(unique(mstrow)) = accumarray(subs(ikeep), prodret(ikeep))./accumarray(subs(ikeep), retSP500(ikeep).^2);
end

end
% %% Analyze flagged
% % Good example: TESTC
% f    = 147;
% d    = '.\data\TAQ';
% dd   = dir(fullfile(d,'*.mat'));
% s    = load(fullfile(d,dd(f).name));
% meta = s.mst{2}(43003,:);
% idx  = 18558662:18559223;
% plot(s.data{2}(idx))
% 
% load('C:\HFbetas\results\flagged.mat')
% baddays = cellfun(@find, baddays,'un',0);
% 
% f = 147;
% disp(baddays{f})
% d       = '.\data\TAQ';
% dd      = dir(fullfile(d,'*.mat'));
% s = load(fullfile(d,dd(f).name));
% scrollto('s.mst{2}',42995)
% idx = 18555526:18555737;
% s.data{2}(idx)
% dates = datenum(double([zeros(numel(idx),3) s.data{1}(idx,:)]));
% fltout(dates,s.data{2}(idx),20,4)
function inspect(ret)
% f = 1;
% s = load(fullfile(d,dd(f).name));
% ret  = s.data.Price(2:end)./s.data.Price(1:end-1)-1;
% ret(setdiff(s.mst.From-1,[-1,0])) = NaN;

[~,pos] = nanmin(ret);
r = find(s.mst.From <= pos & s.mst.To>= pos);
scrollto('s.mst' ,[r ,1])
idx = s.mst.From(r):s.mst.To(r);
s.data.Price(idx)
s.data.Exchange(idx)
end
%% Check stats for returns
function res = dailystats(s,~)
% STEP 1) Selection
inan = selecttrades(s.data);
   
% Select non nans (use prevprice next time)
prices = s.data.Price(~inan);
% Set to NaN overnight
[n, subs] = histc(find(~inan), [1; s.mst.To+1]);

% Daily median price
fill = NaN(1,'single');
nmst = size(s.mst,1);
Med  = accumarray(subs, prices,[nmst,1], @fast_median, fill);

% Min max return
ret              = prices(2:end)./prices(1:end-1)-1;
nnovernight      = true(size(ret));
idx              = setdiff(cumsum(n(1:end-2)),0);
nnovernight(idx) = false;
ret              = ret(nnovernight);
subs             = subs(nnovernight);
% Collect results
res = [accumarray(subs, ret,[nmst,1],@min,fill), accumarray(subs, ret,[nmst,1],@max,fill), Med, n(1:end-1)-1];
end
%% Check stats for returns
function res = avgtimestep(s,cached)
% Number of observations per day
nobs = double(s.mst.To - s.mst.From + 1);
nmst = size(s.mst,1);
% STEP 1) Selection
inan = selecttrades(s.data);

% STEP 2) Bad prices prices are 1.5x or .65x the daily median (net of selection nans)
[~,iprice] = histc(s.data.Price./RunLength(cached.MedPrice,nobs), [.65,1.51]);
iprice     = iprice ~= 1;
inan       = inan | iprice;

% Calculate average intraday unique time step after selection and bad prices
daySubs   = single(RunLength((1:nmst)',nobs));
times     = unique(daySubs(~inan) + single(hhmmssmat2serial(s.data.Time(~inan,:))));
subs      = fix(times);
n         = accumarray(subs, single(1), [nmst,1],   [], ones('single'));
openTime  = accumarray(subs,     times, [nmst,1], @min, ones('single'));
closeTime = accumarray(subs,     times, [nmst,1], @max, ones('single'));
res       = (closeTime - openTime)./(n-1);
end