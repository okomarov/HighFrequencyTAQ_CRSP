function tf = isprobdate(dates)
% ISPROBDATE Checks wether a yyyymmdd date is a problematic one
%
%   NOTE: on those dates trading was half-day, e.g. on the 24th of Dec
probdates = [19931126;19941125;19950703;19951124;19960108;19960705;19961129;19961224;19970703;19971128;19971224;19971226;19981127;19981224;19991126;19991231;20000703;20001124;20010703;20011123;20011224;20020705;20020911;20021129;20021224;20030703;20031128;20031224;20031226;20041126;20051125;20060703;20061124];
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