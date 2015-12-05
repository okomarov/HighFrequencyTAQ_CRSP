function tf = isprobdate(dates)
% ISPROBDATE Checks wether a yyyymmdd date is a problematic one
%
%   NOTE: on those dates trading was half-day, e.g. on the 24th of Dec
probdates = getDates();
try
    tf = ismember(dates,probdates);
catch
    tf = ismember(dates, yyyymmdd2datetime(probdates));
end

% % Incomplete days
% res            = Analyze(testname,[],mst(:, keepflds),[],[],[], badPriceMult);
% res            = loadresults('sampleFirstLast');
% res            = res(res.FirstTime~=0,:);
% [unD, ~, subs] = unique(res.Date);
% idx            = accumarray(subs, res.LastTime,[],@max) < 155900 |...
%                  accumarray(subs, res.FirstTime,[],@min) > 93100;
% incompleteDays = unD(idx);
end

function dates = getDates()
% Closes early at 1 pm on days:
% - following Thanksgiving (4th Thur of November)
% - before/after Independence Day (4th July)
% - Chistmas Eve

years = (1993:year(now))';
cEve  = years*1e4 + 1224;

novDays    = bsxfun(@plus, years*1e4, 1100+(1:30));
isThursday = weekday(yyyymmdd2serial(novDays)) == 5;
daynum     = sum(cumsum(isThursday,2) < 4,2)+2;
afterTGive = years*1e4 + 1100 + daynum;

dates = [...
    19940211 % Snowstorm delays http://blogs.wsj.com/marketbeat/2012/10/29/a-timeline-of-previous-market-shutdowns/
    19950703
    19951218 % Computer system problems http://www.usatoday.com/story/money/2015/07/08/nyse-trading-halted-historical-trading/29866625/
    19960108 % Snowstorm delays http://blogs.wsj.com/marketbeat/2012/10/29/a-timeline-of-previous-market-shutdowns/
    19960705
    19970703
    19971226
    19981026 % Computer system problems http://www.usatoday.com/story/money/2015/07/08/nyse-trading-halted-historical-trading/29866625/
    19991231
    20000703
    20010608 % Computer system problems http://www.usatoday.com/story/money/2015/07/08/nyse-trading-halted-historical-trading/29866625/
    20010703
    20020705
    20020911 % 9/11 memorial http://blogs.wsj.com/marketbeat/2012/10/29/a-timeline-of-previous-market-shutdowns/
    20030703
    20030815 % Power outage https://www.sec.gov/news/testimony/ts102003sec.htm
    20031226
    20060703
    20070703
    20080703
    ];
dates = [dates; cEve; afterTGive];
end