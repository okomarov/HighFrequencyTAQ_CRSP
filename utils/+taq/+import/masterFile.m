function masterFile(fname, outdir, opt)

PREFIX_MASTER = 's_';

% Default case with CSVs
if nargin < 3
    opt.UseTextscan = false;
    opt.ImportFmt   = ['%s %s %s ',repmat('%u8 ',1,9),'%s %f %f %s %u8 %f %*[^\n]'];
    opt.ImportOther = {'Delimiter',',','CommentStyle',{'"','"'}};
end

tb = importTable(fname,opt);

recordsAdd(tb,outdir,PREFIX_MASTER);
end

% Import master file into a table
function tb = importTable(fname,opt)

if opt.UseTextscan
    tb = importTableText(fname,opt);
else
    tb = importTableCsv(fname,opt);
end

% Fill with empty to avoid cat errors
if isempty(tb.DENOM)
    tb.DENOM = repmat(' ', size(tb,1),1);
end
iempty           = cellfun('isempty',tb.NAME);
tb.NAME(iempty)  = {' '};
iempty           = cellfun('isempty',tb.CUSIP);
tb.CUSIP(iempty) = {' '};
end

function tb = importTableText(fname,opt)
% Textscan
fid   = fopen(fname);
clean = onCleanup(@() fclose(fid));
tb    = textscan(fid, opt.ImportFmt,opt.ImportOther{:});

% Conversions
tb(1:3)     = cellfun(@(x) cellstr(x), tb(1:3),'un',0);
tb(4:12)    = cellfun(@(x) logical(x-'0'),tb(4:12),'un',0);
tb{14}      = str2num(tb{14});
tb([15,18]) = cellfun(@(x) uint32(str2num(x)),tb([15,18]),'un',0);
tb{17}      = uint8(tb{17}-'0');

tb = table(tb{:},'VariableNames', opt.VarNames);
end

function tb = importTableCsv(fname,opt)
tb = readtable(fname, 'Format',opt.ImportFmt, opt.ImportOther{:});

% Eventually rename DATEF to FDATE (in some files it changes)
tb.Properties.VariableNames = regexprep(tb.Properties.VariableNames,'(?i)datef','FDATE');

% Conversions
tb = convertColumn(tb, 'logical', {'ETN','ETA','ETB','ETP','ETX','ETT','ETO','ETW','ITS'});
tb = convertColumn(tb, 'int8', 'TYPE');
tb = convertColumn(tb, 'uint32', {'FDATE','UOT'});
tb = convertColumn(tb, 'char', {'ICODE','DENOM'});
end

% Add master records to symbol-specific master files
function tb = recordsAdd(tb, outdir, PREFIX)

% Group records by symbol
[symb,~,subs] = unique(tb.SYMBOL);
tb            = cache2cell(tb,subs);

for ii = 1:numel(symb)
    fname = fullfile(outdir,sprintf([PREFIX '%s'],symb{ii}));

    % Add new records with new dates
    try
        mst = loadMst(fname);
    catch
        saveMst(fname,tb{ii})
        continue
    end

    iold = ismember(tb{ii}.FDATE, mst.FDATE);
    if all(iold)
        continue
    end
    newMst = [mst; tb{ii}(~iold,:)];

    % Sort by CUSIP/DEN and FDATE
    [~,~,newMst.Id] = unique(newMst(:,{'CUSIP','DENOM'}));
    [~,isort]       = sort(uint64(newMst.Id*1e8) + uint64(newMst.FDATE));
    newMst          = newMst(isort,:);

    % Keep unique records with earliest date
    idx       = isfeatchange(newMst(:,[end,2,4:end-1]),[1,3:11,13,14,16,17]);
    newMst    = newMst(idx,:);
    newMst.Id = [];

    if ~isequal(newMst, mst)
        saveMst(fname,newMst)
    end
end
end

function mst = loadMst(fname)
s          = load(fname,'-mat');
mst        = s.(char(fieldnames(s)));
mst.SYMBOL = cellstr(mst.SYMBOL);
mst.NAME   = cellstr(mst.NAME);
mst.CUSIP  = cellstr(mst.CUSIP);
end

function saveMst(fname,mst)
mst.SYMBOL = char(mst.SYMBOL);
mst.NAME   = char(mst.NAME);
mst.CUSIP  = char(mst.CUSIP);
save(fname,'mst','-mat','-v6')
end
