function mapSpconst2mst

% Load spconst
try
    spconst = loadresults('spconst');
catch
    spconst = importDsp500list('.\data\CRSP\');
end

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

% Sample at refdates
refdates = unique(permno2id.Date);
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

% Map to mst
permno2id.Issp500 = ismembIdDate(permno2id.Permno, permno2id.Date, spconst.Permno, spconst.Date);

% Save
issp500 = permno2id;
save(fullfile('.\results',sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'issp500')), 'issp500')

end