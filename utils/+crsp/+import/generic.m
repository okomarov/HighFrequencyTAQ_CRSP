function data = generic(path2zip, outname, fmt, opts)
% GENERIC Import generic CRSP dataset into a table and save as .mat

p              = inputParser();
p.StructExpand = true;
addParameter(p,'VariableDescriptions',[]);
p.parse(opts)

csvfile = unzip(path2zip,tempdir());

fid     = fopen(char(csvfile));
cleanup = onCleanup(@()cleanupFile(fid));

headers = textscan(fid, '%s',1,'Delimiter','');
headers = upperfirst(regexp(headers{1}{1},',','split'));

data = textscan(fid, fmt, 'Delimiter',',');
data = table(data{:}, 'VariableNames', headers);

if ~isempty(p.Results.VariableDescriptions)
    data.Properties.VariableDescriptions = p.Results.VariableDescriptions;
end

filename = sprintf('%s_%s.mat',outname, datestr(now,'yyyymmdd_HHMM'));
save(filename, 'data')
end

function cleanupFile(fid)
fname = fopen(fid);
fclose(fid);
delete(fname)
end
