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
fid      = fopen(char(txtfile));
tline    = fgetl(fid);
isHeader = false;
Desc     = '';
while ischar(tline)
    % Are we parsing the header line? 
    if isHeader
        vnames = textscan(tline, '%s','MultipleDelimsAsOne',true);
        vnames = vnames{1};
        vnames = strrep(vnames,'-','Minus');
        break
    end
    % Empty line, next is header
    if isempty(tline)
        isHeader = true;
    else
        Desc = [Desc, tline];
    end
    tline = fgetl(fid);
end
% Parse data
fmt = ['%u' repmat('%f',1,numel(vnames))];
txt = textscan(fid, fmt,'MultipleDelimsAsOne',true);
fclose(fid);

% Convert percentage returns to decimal
txt(2:end) = cellfun(@(x) x/100, txt(2:end),'un',0);

% Convert into table
data                        = table(txt{:},'VariableNames',['Date'; vnames(:)]);
data.Properties.Description = Desc;

% Delete temporary .txt
delete(txtfile{:})

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