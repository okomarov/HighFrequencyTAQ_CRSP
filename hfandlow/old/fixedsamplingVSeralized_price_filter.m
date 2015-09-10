%% fixedsampling vs sheppard 
% Sheppard does not add any tolerance but uses the <= condition (better)
path2data = '.\data\TAQ';
master = load(fullfile(path2data, 'master'), '-mat');
step = 5/(60*24);
grid = (9.5/24:step:16/24)';

spy = getData(master, 'SPY',19930101,19940101,{'Price'},path2data);
[spydates, idx, subs] = unique(double(spy.Datetime));
spyprices = accumarray(subs, spy.Price,[],@fast_median);


% My sampling
[ps,~,ad] = fixedsampling(spydates, spyprices, grid);
res1 = [ad,ps];

% [ps, ~, ad] = fixedsampling(spydates(idx),spyprices(idx),grid);
% res1 = [ad,ps];
% res2 = [ad2, ps2];

% Sheppards
addpath .\utils\MFE\
time    = rem(spydates,1);
date    = fix(spydates);
unDates = unique(date);
res2 = cell(numel(unDates),1);

for ii = 1:numel(unDates)
    idx = date == unDates(ii);
    [ps2,~,ad2] = realized_price_filter(spyprices(idx),time(idx),...
                                      'unit','fixed',grid);
    res2{ii} = [ad2, ps2];
end
res2 = cat(1,res2{:});

inan = isnan(res1(:,2));
res2 = res2(~inan,:);
res1 = res1(~inan,:);
find(abs(res1(:,2)-res2(:,2)) >= eps) % There are differences


%% chasing the difference

% Verify sheppard's sampling is (a,b]
grid   = 0.1:.1:1;
times  = [0.001, 0.002, 0.09, 0.12, 0.2];
prices = 1:numel(times);

out = zeros(numel(grid),3);
[out(:,3),out(:,1),out(:,2)] = realized_price_filter(prices,times,'unit','fixed',grid);

% ans =
%     0.1    0.09    3
%     0.2    0.20    5
%     0.3    0.20    5
%     0.4    0.20    5
%     0.5    0.20    5
%     0.6    0.20    5
%     0.7    0.20    5
%     0.8    0.20    5
%     0.9    0.20    5
%     1      0.20    5
    
% Test case
dates = [727979.524432870
         727979.538194445];
     
step  = 5/(60*24); % 5 minute step
edges = 727979 + 9.5/24 + step*(37:41); % 40th to 43rd step after 9:30 am
tol   = 1/(60*60*24*1000)-eps;

% Check fixedsampling
nn = histc(dates,edges+tol);

% Check sheppard
prices = 1:numel(dates);
out = zeros(numel(edges),3);
[out(:,3),out(:,1),out(:,2)] = realized_price_filter(prices',dates,'unit','fixed',edges');

isequal(dates(2), edges(end)) % it all lies here! Adding tol within the millisecond changes results.

% I should probably implement a mex file to sample.