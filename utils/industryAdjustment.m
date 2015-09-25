function [signals,ret_monthly] = industryAdjustment(signal, permnos, dates, industry)
% A. "ex-ante industry adjustment" for signals (a-la Novy-Marx (2013) p. 13)
% --------------------------------------------------------------------------


% Intersect TS
idx           = ismember(yyyymmdd2serial(industry.Date), dts_monthly);
industry.Data = industry.Data(idx,:);
industry.Date = industry.Date(idx,:); 

% Match industry XS to ret_monthly
[idx,pos]            = ismember(industry.Permno,constit);
ind_data             = zeros(nobs, nser,'like',industry.Data);
ind_data(:,pos(idx)) = industry.Data(:,idx);

% NaN out XS without industry data
inan                = ~ismember(constit,industry.Permno);
ret_monthly(:,inan) = NaN;
for ii = 1:nsignals
    signals{ii}(:,inan) = NaN;
end

% Lag
nlag     = 1;
ind_data = [zeros(nlag,nser); ind_data(1:end-nlag,:)];

% Number of industries and subscripts
[unind,~,avgsubs] = unique(ind_data);
num_industries    = numel(unind);

% Loop for each signal
for hh = 1:nsignals
    signal  = signals{hh};
    nobs    = size(signal,1);
    % TS of signals averaged by industry
    ind_avg = NaN(nobs,num_industries);
    for ii = 1:nobs
        for jj = 1:num_industries
            idx            = ind_data(ii,:) == unind(jj);
            ind_avg(ii,jj) = nanmean(signal(ii,idx)); % [Do!] Perhaps in future code VW vs EW choice to be consistent with b part
        end
    end
    % Set to NaN class-less industry
    ind_avg(:,unind == 0) = NaN;
    
    % Demean signal by its industry average
    rowsubs = repmat((1:nobs)',1,nser);
    iavg    = sub2ind([nobs,nser], rowsubs(:), avgsubs);
    ind_avg = reshape(ind_avg(iavg), [nobs,nser]);
    signal  = signal - ind_avg;
    
    signals{hh} = signal;
end