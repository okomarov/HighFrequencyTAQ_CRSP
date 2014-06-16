function res = Analyze(fun, varnames, cached, debug)

% ANALYZE Executes specified fun in parallel on the whole database (all .mat files) 
% 
%   ANALYZE(FUN, VARNAMES) FUN should a string with the name of one of 
%                          the following sub-functions:
%                               - 'dailystats'            
%                               - 'badprices'
%                               - 'avgtimestep'
%                          VARNAMES is a cell-array of strings (or string)
%                          with the VarNames of the dataset with the results
%
%   ANALYZE(..., CACHED) Some FUN might require pre-cached results which where
%                        run on the whole database.
%                        Check the specific sub-function for the format of the 
%                        needed CACHED results.
%   ANALYZE(..., DEBUG) Run execution sequentially, i.e. not in parallel, to be 
%                       able to step through the code in debug mode.


if nargin < 3; cached = [];    end
if nargin < 4; debug  = false; end   

% Simply call the specific subroutine
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
    dd  = dir(fullfile(d,'T*.mat'));
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
        res{f} = fhandle(s, [cached(f), f]);
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
%% Check stats for returns
function res = dailystats(s,~)
% DAILYSTATS Relevant for median calculation

% STEP 1) Selection
inan = selecttrades(s.data);
  
% Select non nans
prices = s.data.Price(~inan);
% Set to NaN overnight
[n, subs] = histc(find(~inan), [1; s.mst.To+1]);

% Daily median price
fill = NaN('single');
nmst = size(s.mst,1);
Med  = accumarray(subs, prices,[nmst,1], @fast_median, fill);

% Min max return
ret              = prices(2:end)./prices(1:end-1)-1;
nnovernight      = true(size(ret));
idx              = setdiff(cumsum(n(1:end-2)),0);
nnovernight(idx) = false;
ret              = ret(nnovernight);
subs             = subs(nnovernight);
% Collect results (note, nrets can be negative if it was only one price, which was selected out)
res = [accumarray(subs, ret,[nmst,1],@min,fill), accumarray(subs, ret,[nmst,1],@max,fill), Med, n(1:end-1)-1];
end
%% Count bad prices
function res = badprices(s, cached)
cached = cached{1};
% Flag bad days if number of intraday bad prices is > dailycut
dailycut = .5;

% Number of observations per day
nobs = double(s.mst.To - s.mst.From + 1);

% STEP 1) Selection
inan = selecttrades(s.data);

% STEP 2) Bad prices are < than .65x daily median or > than 1.5x daily median
[~,ibadprice] = histc(s.data.Price./RunLength(cached.MedPrice,nobs), [.65,1.51]);
ibadprice     = ibadprice ~= 1;

% STEP 3) Bad days
res = accumarray(RunLength((1:size(s.mst,1))',nobs), inan | ibadprice) > ceil(dailycut*nobs);
end
%% Check stats for returns
function res = avgtimestep(s,cached)
cached = cached{1};

% Number of observations per day
nobs = double(s.mst.To - s.mst.From + 1);
nmst = size(s.mst,1);
% STEP 1) Selection
inan = selecttrades(s.data);

% STEP 2) Bad prices are < than .65x daily median or > than 1.5x daily median
[~,ibadprice] = histc(s.data.Price./RunLength(cached.MedPrice,nobs), [.65,1.51]);
ibadprice     = ibadprice ~= 1;
inan          = inan | ibadprice;

% STEP 3) Take median price if same timestamp
mstrow = RunLength((1:nmst)',nobs);
mstrow = mstrow(~inan);
times  = single(unique(mstrow + hhmmssmat2serial(s.data.Time(~inan,:))));

% Calculate average intraday unique time step after selection, bad prices and consolidation
subs      = uint32(fix(times));
l         = ones('single');
n         = accumarray(subs,      l, [nmst,1],   [], l);
openTime  = accumarray(subs,  times, [nmst,1], @min, l);
closeTime = accumarray(subs,  times, [nmst,1], @max, l);
res       = (closeTime - openTime)./(n-1);
end
%% Check how many bad prices
function res = cleanprices(s,cached)
cached = cached{1};

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
%% Sampling
function res = sample(s,cached)
nfile  = cached{end};
cached = cached{1};
% Number of observations per day
nobs   = double(s.mst.To - s.mst.From + 1);
nmst   = size(s.mst,1);

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

% STEP 5) Clean prices - carried out in the median consolidation
% s.data.Price(inan) = NaN;

if ~all(inan)

    % STEP 6) Median prices for same timestamps
    mstrow           = RunLength((1:nmst)',nobs);
    [unTimes,~,subs] = unique(mstrow(~inan) + hhmmssmat2serial(s.data.Time(~inan,:)));
    price            = accumarray(subs, s.data.Price(~inan),[],@fast_median);
        
    % STEP 7) Sample on fixed grid (easier to match sp500) 
    grid  = (9.5/24:5/(60*24):16/24)';
    ngrid = numel(grid);
    price = fixedsampling(double([unTimes, price]), 'Previous', grid);
    idx   = fix(price(:,1));
    dates = yyyymmdd2serial(double(s.mst.Date(idx))) + rem(price(:,1),1);
    
    % Save
    res   = [];
    fname = fullfile('.\data\TAQ\sampled',sprintf('S5m_%04d.mat',nfile));
    data  = dataset({dates,'Datetime'},{price(:,2),'Price'});
    
    imst     = RunLength(idx);
    mst      = [s.mst(imst,'Id'), cached(imst, 'UnID'), s.mst(imst,'Date')];
    mst.From = (1:ngrid:size(data,1))';
    mst.To   = (ngrid:ngrid:size(data,1))';
    
    ids = s.ids;
    save(fname, 'data','ids', 'mst')
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