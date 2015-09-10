function res = mapUnid2mst(mst, symbols)
% MAPUNID2MST Maps UnIDs to master records with Id and Date and Symbol

if isa(mst,'dataset'), mst = dataset2table(mst); end

% Replace p with PR and remove . from TAQ symbols
symbols = regexprep(symbols,'p','PR');
symbols = regexprep(symbols,'\.','');

% Get taq2crsp ready
taq2crsp = loadresults('taq2crsp');
taq2crsp = sortrows(taq2crsp,{'symbol','datef'});

% Keep only cusip match (score 10)
taq2crsp = taq2crsp(taq2crsp.score == 10,:);

% % Remove empty cusips
% taqmaster = loadresults('TAQmaster');
% taqmaster = taqmaster(~cellfun(@isempty, taqmaster.CUSIP),:);
% taqmaster = sortrows(taqmaster,{'CUSIP','FDATE'});
% idx = isfeatchange(taqmaster(:,{'CUSIP','SYMBOL','FDATE'}),3);
% taqmaster = taqmaster(idx,:);

% Preallocate
unID = zeros(size(mst,1),1,'uint16');

% LOOP for each symbol in mst
for ii = 1:numel(symbols)
    symbol   = symbols{ii};
    % Extract records from taq2crsp corresponding to TAQ's symbol
    isymbol  = strcmpi(symbol,taq2crsp.symbol);
    tmp      = taq2crsp(isymbol, {'ID','datef'});
    if isempty(tmp),fprintf('%d\n',ii),continue,end
    imst     = find(mst.Id == ii);
    % Find to which date interval the records belong
    [~,itmp] = histc(mst.Date(imst), [tmp.datef; 99999999]);
    nnzero   = itmp ~= 0;
    % Assgin unique ID to mst
    unID(imst(nnzero)) = tmp.ID(itmp(nnzero));
end
% Unmatched
unID(unID == 0) = intmax('uint16');
res = [mst table(unID,'VarNames','UnID')];
save(fullfile('.\results\', sprintf('%s_%s.mat', datestr(now,'yyyymmdd_HHMM'),'uniqueID')), 'res')
end