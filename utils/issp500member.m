function tf = issp500member(tb)
% ISSP500MEMBER Checks which Id - Date pair are sp500 members
%
%   ISSP500MEMBER(TB) TB is a table with Id and yyyymmdd Date 


if isa(tb,'dataset')
    tb = dataset2table(tb);
end

% Load sp500 membership table
spconst = loadresults('spconst');

% Sample at refdates
refdates = unique(tb.Date);
spconst  = sampledates(spconst.Panel, refdates,false);

% Stack back into tall table
vnames  = getVariableNames(spconst);
spconst = stack(spconst, vnames(2:end), ...
                'IndexVariableName', 'Permno',...
                'NewDataVariableName','Issp500');

% Keep only SP500 members
spconst = spconst(spconst.Issp500 ~= 0,:);

% Convert literal xPermno to numeric id
Permno         = char(spconst.Permno);
Permno         = Permno(:,2:end)';
Permno         = textscan([Permno; repmat(' ',1,size(Permno,2))],'%u32');
spconst.Permno = Permno{1};

% Load Permno to Id mapping
try
    permno2id = loadresults('masterPermnoSP500');
catch
    % Load Permno - date pairs
    permno2id = loadresults('masterPermno');
   
    % Filter out non members
    idx       = ismember(permno2id.Permno, spconst.Id);
    permno2id = permno2id(idx,:);
    
    % Cache results
    fname = fullfile('.\results\',sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'masterPermnoSP500'));
    save(fname,'id2permno')
end

% Map membership to Id
idx       = ismembIdDate(permno2id.Permno, permno2id.Date,spconst.Permno, spconst.Date);
permno2id = permno2id(idx,:);

% Filter out tb
tf = ismembIdDate(tb.Id, tb.Date, permno2id.Id, permno2id.Date);
end