function avg = retbyday(dt,ret, freq)
% RETBYDAY Average return by day within a specifed grouping frequency
%
%   dt should be yyyymmdd

if nargin < 3
    freq = 'month';
end
cases = {'month','week'};
icase = strncmpi(freq, cases, numel(freq));
freq  = cases{icase};

switch freq
    case 'month'
        [~,~,idate] = unique(dt/100);
        ndays       = 20;

        [num, n] = deal(zeros(ndays,size(ret,2)));
        allnan   = true(ndays,size(ret,2));

        for ii = 1:max(idate)
            idx      = idate == ii;
            retslice = ret(idx,:);
            inan     = isnan(retslice);

            retslice = nan2zero(retslice([1:10,end-9:end],:));
            inan     = inan([1:10,end-9:end],:);

            % Mean as S(sum)/S(n)
            num    = num + retslice;
            n      = n + ~inan;
            allnan = allnan + inan;
        end
        avg         = num./n;
        avg(allnan) = NaN;

    case 'week'
        ndays    = 5;
        avg = zeros(ndays,size(ret,2));
        idate = weekday(yyyymmdd2serial(dt))-1;

        for ii = 1:max(idate)
            idx      = idate == ii;
            retslice = ret(idx,:);
            avg(ii,:) = nanmean(retslice,1);
        end

    otherwise
        error('error:unrecognizedFreq','FREQ not recognized.')

end
end
