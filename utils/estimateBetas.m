function betas = estimateBetas(lookback, freq, useon, useproxy)
if nargin < 1 || isempty(lookback), lookback = 1;     end
if nargin < 2 || isempty(freq),     freq     = 5;     end
if nargin < 3 || isempty(useon),    useon    = true; end
if nargin < 4 || isempty(useproxy), useproxy = false; end

writeto = '.\results\';

try
    fprintf('%s: loading betacomponents at %d min.\n', mfilename, freq)
    name  = matname('betacomponents',freq, useon, useproxy);
    betas = loadresults(name);
catch
    fprintf('%s: betacomponents not found. Estimating...\n', mfilename)
    
    % Sample if data doesn't exist
    path2data = sprintf('.\\data\\TAQ\\sampled\\%dmin', freq);
    if exist(path2data,'dir') ~= 7 || numel(dir(path2data)) <= 2
        fprintf('%s: sampling data at %d min.\n', mfilename, freq)
        step    = freq/(60*24);
        grid    = (9.5/24:step:16/24)';
        fmtname = sprintf('S%dm_%%04d.mat',freq);
        sampleData(grid, path2data, fmtname);
    end 
        
    % Use self-built sp500 proxy
    if useproxy
        name = sprintf('sp500proxy%dm',freq);
        try
            sp500 = loadresults(name);
        catch
            fprintf('%s: creating ssp500proxy at %d min.\n', mfilename, freq)
            sp500 = sp500intraday(path2data);
        end
    else
        % Use spyders
        sp500 = getSpy(freq);
    end
    
    % SP500 ret 
    spret = [sp500.Datetime, [NaN; diff(log(sp500.Price))]];
    ion   = [false; diff(rem(sp500.Datetime,1)) >= 0];
    
    % Add overnight return or discard
    if useon
        fprintf('%s: adding overnight returns to the index.\n', mfilename)
        ion               = find(~ion);
        reton             = loadresults('return_intraday_overnight');
        spreton           = reton(reton.Permno == 84398,:);
        [idx,pos]         = ismember(serial2yyyymmdd(spret(ion,1)),spreton.Date);
        spret(ion(idx),2) = spreton.RetCO(pos(idx));
    else
        spret = spret(ion,:);
    end
        
    % Cache SP500 returns by days
    fprintf('%s: caching index returns by days.\n', mfilename)
    load(fullfile(path2data,'master'),'-mat','mst');
    nfiles = max(mst.File);
    cached = cell(nfiles,2);
    
    dates           = fix(spret(:,1));
    [spdays,~,subs] = unique(dates,'stable');
    spret           = mat2cell(spret(:,2), accumarray(subs,1),1);
    
    unMstDates = accumarray(mst.File, mst.Date,[],@(x){yyyymmdd2serial(unique(x))});
    for ii = 1:nfiles
        [~,pos]      = ismember(unMstDates{ii}, spdays);
        nnzero       = pos ~= 0;
        isp          = ismember(spdays, unMstDates{ii});
        cached(ii,:) = {spret(isp) spdays(pos(nnzero))};
    end
    
    % Cache overnight returns
    % Note: the match is on Date -Id rather than Date - Permno to avoid the
    %       duplication of overnight returns that comes from the 
    %       Date - Permno mapping to Date - Permno/Id
    if useon
        fprintf('%s: adding overnight returns to all other series.\n', mfilename)
        reton.File      = zeros(size(reton,1),1,'uint16');
        [idx,pos]       = ismembIdDate(reton.Id,reton.Date, mst.Id, mst.Date);
        reton.File(idx) = mst.File(pos(idx));
        reton           = reton(reton.File ~= 0,{'Permno','Date','RetCO','File'});
        cached          = [cached,...
                           accumarray(reton.File,(1:size(reton))',[],@(x) {reton(x,{'Permno','Date','RetCO'})})];
    end
    
    % Calculate beta components: sum(r*benchr) and sum(benchr^2)
    fprintf('%s: creating betacomponents at %d min.\n', mfilename, freq)
    [betas,filename] = Analyze('betacomponents', [], cached, path2data);
    
    % Rename to append the sampling frequency
    name        = regexp(filename,'\w+?(?=\.mat)','match','once');
    name        = [matname(name,freq, useon, useproxy),'.mat'];
    newfullname = fullfile(writeto, name);
    movefile(fullfile(writeto,filename), newfullname);
end
fprintf('%s: calculating betas with %d day lookback.\n', mfilename, lookback)

% Sort (composite key for speed) and create subs
key        = uint64(betas.Permno) * 1e8 + uint64(betas.Date);
[~,isrt]   = sort(key);
betas      = betas(isrt,:);
[~,~,subs] = unique(betas.Permno);

% Beta
num        = accumarray(subs, betas.Num, [], @(x) runsum(lookback,x));
den        = accumarray(subs, betas.Den, [], @(x) runsum(lookback,x));
betas.Beta = cat(1,num{:})./cat(1,den{:});
end

function name = matname(name, freq, useon, useproxy)
if useon,    useon    = 'on';   else useon    = ''; end
if useproxy, useproxy = 'prx';  else useproxy = ''; end
name = sprintf('%s%dm%s%s', name, freq, useon, useproxy);
end