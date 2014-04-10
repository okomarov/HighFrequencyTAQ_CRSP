%% Common vars
d     = '.\data\SP500';
start = datenum('31/12/1992','dd/mm/yyyy');
opts  = {'Delimiter',',','ReadVarNames',1,'CommentStyle',{'"','"'}};

%% Load: SP500 series for dates vector and members count cross check
% Note: there are holiday effects in the member counts, i.e. sometime some constituent ceases to be a SP
% member before a holiday and the new member is added after the holiday. Thus, to get proper counts you
% need to refer to the proper vector of dates
SP500crsp = dataset('File',fullfile(d,'dsp500.csv'),opts{1:4});
SP500crsp = replacedata(SP500crsp,@yyyymmdd2serial,'caldt');

% Keep those with enddate > 31/12/1992
SP500crsp = SP500crsp(SP500crsp.caldt > start,:);
%% Load: costituents SP500
const = dataset('File',fullfile(d,'dsp500list.csv'),opts{1:4});
const = replacedata(const,@yyyymmdd2serial,{'start','ending'});

% Keep those with enddate > 31/12/1992
const = const(const.ending > start,:);

% Dates for pivoting (get rid of the weekends)
dates = SP500crsp.caldt;

% Pivot it
[idConst,~,iConstp] = unique(const.PERMNO);
nid        = numel(idConst);
constp     = false(numel(dates),nid);
for ii = 1:size(const,1)
    constp(:, iConstp(ii)) = constp(:, iConstp(ii)) | dates >= const.start(ii) & dates <= const.ending(ii);
end

const = struct('constdt',const,'constp',constp,'ids',idConst,'dates',dates);
%% Check: SP500 member daily counts 
% Count number of members
counts = sum(const.constp,2);

% % Alternative count (SAME RESULT as previous)
% counts = zeros(numel(const.dates),1);
% for ii = 1:numel(const.dates)
%     checkdate = const.dates(ii);
%     counts(ii) = numel(unique(const.constdt.PERMNO(checkdate >=  const.constdt.FROM & checkdate <= const.constdt.TO)));
% end

% Check against reference cuonts
isequal(counts, SP500crsp.totcnt)
%% Load: number of shares
% 5 secs to import
nshares = dataset('File',fullfile(d,'dseshares.csv'),'Format','%f%f%f%f%*[^\n]',opts{:});
nshares = replacedata(nshares,@yyyymmdd2serial,{'SHRSDT','SHRENDDT'});

% Keep those with enddate > 31/12/1992
nshares = nshares(nshares.SHRENDDT > start,:);

% Pivot only constituents - 40 secs
tf         = ismember(nshares.PERMNO,idConst);
nsharesred = nshares(tf,:);
[idNshares,~,iNsharesp] = unique(nsharesred.PERMNO);
nid        = numel(idNshares);
nsharesp   = zeros(numel(dates),nid);
for ii = 1:size(nsharesred,1)
    idx = dates >= nsharesred.SHRSDT(ii) & dates <= nsharesred.SHRENDDT(ii);
    nsharesp(idx, iNsharesp(ii)) = nsharesred.SHROUT(ii);
end
nshares = struct('nsharesdt',nshares,'nsharesp',nsharesp,'ids',idNshares,'dates',dates);
%% Check: if member, has shrout? 
% Member but no shrout
pos = find(nshares.nsharesp==0 & const.constp);
% Convert to row/col subs
szp   = size(nshares.nsharesp);
[r,c] = ind2sub(szp, pos);

% Find first position before and after the member without shrout
beforeAfter      = zeros(size(pos,1),2);
dfpos            = diff([pos;prod(szp)]);
beforeAfter(1,:) = [pos(1)-1, pos(find(dfpos > 1,1,'first'))+1];
for ii = 2:numel(pos)
    if pos(ii) == pos(ii-1)+1
        beforeAfter(ii,:) = beforeAfter(ii-1,:);
    else
        beforeAfter(ii,:) = [pos(ii)-1, pos(find(dfpos(ii:end)>1,1,'first')+ii-1) + 1];
    end
end

% Extend (either before or after)
extendval = nshares.nsharesp(beforeAfter);
whichcol  = (extendval == 0) * [1;2];
extendval = max(extendval,[],2);
for ii = 1:numel(pos)
    col = whichcol(ii);
    if col == 2
        idx = pos(ii):beforeAfter(ii,col)-1;
    else
        idx = beforeAfter(ii,col)+1:pos(ii);
    end
    nshares.nsharesp(idx) = extendval(ii);
end
isempty(find(nshares.nsharesp==0 & const.constp,1)) % check
% Inspect
% n = ii;
% scrollto('nshares.nsharesp',pos(n))
% scrollto('const.constp',pos(n))
% scrollto('const.dates',r(n))
%% Load: PERMNO to NCUSIP (historical CUSIP)
fmt        = '%f%f%f%f%s%s%s%s%f%f%f%f%s%f%f%f';
stocknames = dataset('File',fullfile(d,'stocknames.csv'), 'Format',fmt,opts{:});
stocknames = replacedata(stocknames,@yyyymmdd2serial, {'NAMEDT','NAMEENDDT','ST_DATE','END_DATE'});
                                 
% Keep only constituents and NCUSIPs valid after 1992/12/31
stocknames = stocknames(stocknames.NAMEENDDT > start & ismember(stocknames.PERMNO, const.ids),:);

%% Load CSV TAQ master files
load('data')
% List all TAQ .csv master files
d         = '.\data\TAQ\raw';
list      = unzip(fullfile(d,'mast_csv.zip'),d);
list      = sort(list);
nfiles    = numel(list);
TAQmaster = cell(nfiles,1);
fmt       = ['%s %*s %s ',repmat('%*u8 ',1,9),'%*s %f %*f %*s %*u8 %f %*[^\n]'];

% LOOP by file - 20 sec
for ii = 1:nfiles
    disp(list{ii})
    % Read in the whole master file
    TAQmaster{ii} = dataset('File',list{ii}, 'Format',fmt, opts{:});
    % Eventually rename DATEF to FDATE (in some files it changes)
    TAQmaster{ii}.Properties.VarNames = regexprep(TAQmaster{ii}.Properties.VarNames,'(?i)datef','FDATE');
    % Change date to Matlab serial number
    TAQmaster{ii} = replacedata(TAQmaster{ii},@yyyymmdd2serial,'FDATE');
end

% Concatenate everything - 20 sec
TAQmaster = cat(1,TAQmaster{:});

% Keep first (by default from 2013a) changes of the tuple symbol, cusip, shrout - 5 sec
TAQmaster = unique(TAQmaster,{'SYMBOL','CUSIP','SHROUT'});

% Extract the 8-CUSIP
taqcusips = regexp(TAQmaster.CUSIP,'[^0]\w*(?=\d{4}$)','match','once');
    
% Sort
TAQmaster = sortrows(TAQmaster,{'SYMBOL','FDATE'});

% Clean
delete(list{:})

%% Checkpoint
clearvars -except SP500crsp const nshares stocknames opts start TAQmaster taqcusips
save('data.mat')
%% Daily processing (Redo)

% Load TAQ master data
d         = '.\data\TAQ';
mstdata   = load(fullfile(d,'\master'),'-mat');
mstdata   = mstdata.mst;

% List all TAQ .csv master files
list      = dir(fullfile(d,'raw','*.csv'));
nfiles    = numel(list);
TAQmaster = cell(nfiles,1);
fmt       = ['%s%*s%s',repmat('%*u8',1,9),'%*s%f%*f%*s%*u8%f'];
% LOOP by month
for ii = 1:nfiles
    % Read in the master file
    TAQmaster{ii} = dataset('File',fullfile(d,'raw',list(ii).name), 'Format',fmt, opts{:});
    TAQmaster{ii} = replacedata(TAQmaster{ii},@yyyymmdd2serial,'FDATE');
    taqcusips     = regexp(TAQmaster{ii}.CUSIP,'[^0]\w*(?=\d{4}$)','match','once');
    
    % Select active SP500 dates of the year/month of the master file
    [yy,mm] = datevec(TAQmaster{ii}.FDATE(1));
    enddate = datenum(yy,mm+1,0);
    dates   = const.dates(const.dates <= enddate);
    ndates  = numel(dates);
    % LOOP by day
    for jj = 1:ndates
        jj
        % Index constituents for the date
        iconst = const.constp(const.dates == dates(jj),:);
        % Select the active PERMNOs 
        icusip = ismember(stocknames.PERMNO, const.ids(iconst)) &...
                 stocknames.NAMEDT <= dates(jj) & dates(jj) <= stocknames.NAMEENDDT ;
        % Fetch the corresponding NCUSIP
        cusips = stocknames.NCUSIP(icusip);
        ncusips = numel(cusips);
        if ncusips ~= numel(unique(cusips))
            warning('Less unique NCUSIPs than PERMNOs!')
        end
        % Index the TAQ cusips
        [itaq,subs] = ismember(taqcusips, cusips);
        taqsymb = zeros(ncusips,1);
        % Retrieve the last active SYMBOL for the CUSIP
        for zz = 1:numel(cusips)
            taqsymb(zz) = find(TAQmaster{ii}.FDATE <= dates(jj) & subs == zz,1,'first');
        end
        taqsymb = TAQmaster{ii}.SYMBOL(taqsymb);
        if ncusips ~= numel(unique(taqsymb))
            warning('Less unique TAQSYMBOLS than NCUSIPS')
        end
        % Retrieve ID in master data
        [~,taqid] = ismember(taqsymb, mstdata{1});
        % Check which mat files to load
        mstsel = mstdata{2}(ismembc(mstdata{2}(:,1),sort(uint32(taqid))) & mstdata{2}(:,2) == serial2yyyymmdd(dates(jj)),:);
        if numel(unique(mstsel(:,end))) > 1
            error
        end
    end
end
%% (OLD) Load SP500 costituents: GVKEY, IID and from thru
% d           = '.\data\CCM_SP500';
% fid         = fopen(fullfile(d,'sp500constituents.csv'));
% costituents = textscan(fid,'%f%f%*f%s%s%*[^\n]','HeaderLines',1,'Delimiter',',');
% fclose(fid);
% 
% % Index empty 'to' dates
% iEmpty = cellfun('isempty',costituents{4});
% 
% % Convert to dataset
% costituents = dataset({costituents{1},'GVKEY'},{costituents{2},'IID'},...
%                       {datenum(costituents{3},'yyyymmdd'),'from'},...
%                       {datenum(costituents{4},'yyyymmdd'),'to'});
% 
% % Replace indexed empty with inf
% costituents.to(iEmpty) = inf;