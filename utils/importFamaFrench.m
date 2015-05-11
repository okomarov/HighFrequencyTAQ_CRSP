function data = importFamaFrench(zipname, outdir)
% IMPORTFAMAFRENCH Imports datasets from Kenneth French's Data Library webpage
%
%   IMPORTFAMAFRENCH(ZIPNAME) Imports into a table the dataset specified
%                             by ZIPNAME.
%
%   IMPORTFAMAFRENCH(...,OUTDIR) Specify in which folder to save the
%                                imported data.
%                                By default it will be saved under
%                                '.\ZIPNAME.mat'
%
%   IMPORTFAMAFRENCH() Lists available datasets, their ZIPNAMEs and the
%                      description.
%
%
% See <a href="http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html">Fama French data</a>

url = 'http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/';

% List available datasets on Kennet's page
if nargin < 1
    list = listAvailableData(url);
    data = list;
    return
end
if nargin < 2
    outdir = '';
end

% Try to load from web
try
    txtfile = unzip([url, 'ftp/', zipname], tempdir);
catch
    error('Problem with importing')
end

% Parse .txt dataset:
%   Usually has few lines of description followed by one blank line, one
%   line of variable names, and then data. Uses multiple whitespaces as
%   delimiter.
fid   = fopen(char(txtfile));
clean = onCleanup(@() cleanupFcn(fid,txtfile{:}));
tline = fgetl(fid);
Desc  = '';
while ischar(tline)
    if isempty(tline)
        continue
    end
    
    % First line of data?
    tmp = textscan(tline, '%d','MultipleDelimsAsOne',true);
    tmp = tmp{1};
    if ~isempty(tmp)
        % Parse previous line for variable names
        vnames   = textscan(prevline, '%s','MultipleDelimsAsOne',true);
        vnames   = vnames{1};
        ncoldata = numel(tmp);
        % Create variable names
        if ncoldata == numel(vnames)+1 && vnames{end}(end) ~= '.'
            vnames = strrep(vnames,'-','Minus');
        else
            vnames = matlab.lang.makeUniqueStrings(repmat({'Var'},1,ncoldata-1));
        end
        vnames = ['Date'; vnames(:)];
        
        % Rewind to beginning of data
        fseek(fid, prevpos, 'bof');
        break
    else
        % Store description
        Desc = [Desc, tline];
    end
    
    prevpos  = ftell(fid);
    prevline = tline;
    tline    = fgetl(fid);
end

% Parse data
fmt = ['%u' repmat('%f',1, ncoldata-1)];
txt = textscan(fid, fmt,'MultipleDelimsAsOne',true);

% Convert into table
data                        = table(txt{:},'VariableNames', vnames);
data.Properties.Description = Desc;

% Save
[~,name] = fileparts(zipname);
save(fullfile(outdir,[name,'.mat']), 'data')
end

function list = listAvailableData(url)
% Import webpage
str          = webread([url,'/data_library.html']);
% Parse links
[Link,enpos] = regexp(str, '<a href ?= ?"([^"]*)">','tokens','end');
Link         = cat(1,Link{:});
% Get zipnames
idx          = ~cellfun('isempty',regexp(Link,'.zip','once'));
Zipname      = strrep(Link(idx),'ftp/','');
enpos        = enpos(idx);
list         = table(Zipname);
% Extract description
for ii = 1:numel(enpos)
    p                      = enpos(ii);
    tmp                    = regexprep(str(p+1:p+150),'(<[ba]>|</b>)','');
    list.Description(ii,1) = strtrim(regexp(tmp, '.*(?=</a>)','match'));
end
end

function cleanupFcn(fid,fname)
% Cleanup performed at end or error
fclose(fid);
delete(fname);
end