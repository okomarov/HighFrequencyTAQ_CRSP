function tf = issp500member(tb)
% ISSP500MEMBER Checks which UnID - Date pair are sp500 members
%
%   ISSP500MEMBER(TB) TB is a table with UnID and yyyymmdd Date 


if isa(tb,'dataset')
    tb = dataset2table(tb);
end
% Load sp500 membership table
spconst = loadresults('spconst');

try
    taq2crsp = loadresults('taq2crsp_sp500');
catch
    % Load unid2permno
    taq2crsp = loadresults('taq2crsp');
   
    % Keep only match score 10
    taq2crsp = taq2crsp(taq2crsp.score == 10,:);
        
    % Filter out non members
    taq2crsp   = dataset2table(taq2crsp(ismember(taq2crsp.permno, spconst.Id),:));
    
    % Drop tickers with suffixes
    [~,~,subs] = unique(taq2crsp.permno);
    N          = numel(subs);
    idrop      = false(N,1);
    
    % LOOP for each permno
    for ii = 1:N
        % Index permno and corresponding tickers
        isubs   = subs == ii;
        symbols = taq2crsp.symbol(isubs);
        nsymb   = nnz(isubs);
        
        if nsymb > 1
            idx = false(nsymb, 1);
            for jj = 1:nsymb
                if all(idx), break,end
                % Mark to drop if substring
                idx = idx | ~cellfun('isempty', regexp(symbols, sprintf('^%s\\w+',symbols{jj}),'once'));
            end
            if any(idx), idrop(isubs) = idx; end
        end
    end
    taq2crsp = taq2crsp(~idrop,:);

    % Cache results
    fname = fullfile('..\results\',sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'taq2crsp_sp500'));
    save(fname, 'taq2crsp')
end

% Map UnID to permno
[idx,pos]      = ismember(tb.UnID, taq2crsp.ID);
tbsmall        = tb(idx,'Date');
tbsmall.Permno = taq2crsp.permno(pos(idx));

% Remove overlapping
[un,~,subs] = unique(tbsmall);
overlap     = un(accumarray(subs,1) > 1,:);
% Speed up with composite key
keyA        = uint64(tbsmall.Permno) * 1e8 + uint64(tbsmall.Date);
keyB        = uint64(overlap.Permno) * 1e8 + uint64(overlap.Date);
ioverlap    = ismember(keyA, keyB);
tbsmall     = tbsmall(~ioverlap,:);
idx(idx)    = idx(idx) & ~ioverlap;

% taq2crsp(ismember(taq2crsp.permno, unique(overlap.Permno)),:)

% Pivot tb
warning off MATLAB:table:ModifiedVarnames
nrows       = size(tbsmall,1);
tbsmall.Val = ones(nrows,1,'uint8');
tbpanel = unstack(tbsmall,'Val','Permno');
permnos = unique(tbsmall.Permno);
warning on MATLAB:table:ModifiedVarnames

% Intersect permnos
tbnames = getVariableNames(tbpanel); 
spnames = getVariableNames(spconst.Panel); 
[~,itb,isp] = intersect(tbnames, spnames);
tbpanel = table2array(tbpanel(:,itb));
sppanel = table2array(spconst.Panel(:,isp));

% Intersect/expand dates
refdates       = tbpanel(:,1);
alldates       = union(sppanel(:,1), refdates);
[~,pdates]     = ismember(sppanel(:,1),alldates);
data           = NaN(numel(alldates),size(sppanel,2)-1);
data(pdates,:) = sppanel(:,2:end);
% Fill in previous value and get reference dates only
data           = nanfillts(data,1);
sppanel        = [refdates data(ismember(alldates, refdates),:)];

tbpanel(~sppanel) = 0;
tbpanel = array2table(tbpanel,'VariableNames',tbnames);
tbpanel = stack(tbpanel,tbnames(2:end),'NewDataVariableName','Val','IndexVariableName','Permno');
tbpanel = tbpanel(tbpanel.Val == 1,1:2);
[~,pos] = ismember(tbpanel.Permno,tbnames(2:end));
tbpanel.Permno = permnos(pos);

% Speed up with composite key
keyA   = uint64(tbsmall.Permno) * 1e8 + uint64(tbsmall.Date);
keyB   = uint64(tbpanel.Permno) * 1e8 + uint64(tbpanel.Date);
ismall = ismember(keyA, keyB);

tf      = false(size(tb,1),1);
tf(idx) = ismall;

end