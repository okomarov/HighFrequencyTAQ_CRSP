function out = formatResults(coeff, se, pval)

% Get optional arguments
fmtCoeff = '%.2f';
fmtSe    = '\\tiny(%.3f)';

% Pre-allocate
sz  = size(coeff);
out = cell(sz(1)*2, sz(2));

% Format coefficients and se
cellcoeff = arrayfun(@(x)sprintf(fmtCoeff,x), coeff, 'un', 0);
out(2:2:end,:) = arrayfun(@(x)sprintf(fmtSe,x), se, 'un', 0);

% Add *, **, *** for significance at 10, 5 and 1%
a = [-inf, 0.01, 0.05, 0.1];
for ii = 1:3
    idx = a(ii) < pval & pval <= a(ii+1);
    repstr = sprintf('$0\\\\textsuperscript{%s}',repmat('\*',1,4-ii));
    cellcoeff(idx) =  regexprep(cellcoeff(idx), '.*',repstr); 
end
out(1:2:end,:) = cellcoeff;

% Clean out NaNs
inan = logical(kron(isnan(coeff),ones(2,1)));
out(inan) = {[]};

% Add fixed width separation
len    = cellfun(@length, out);
maxlen = max(len(:));
N      = numel(out);
for ii = 1:N-sz(1)*2
    out{ii} = sprintf('%s%s &',repmat(' ',1, maxlen-len(ii)), out{ii});
end

% Add endline
for ii = N:-1:N-sz(1)*2+1
    out{ii} = sprintf('%s%s \\\\',repmat(' ',1, maxlen-len(ii)), out{ii});
end

% Add column headers
cHeaders = {'HF','LF','HF-LF'};
ncols = sz(2)/numel(cHeaders);
for ii = 1:numel(cHeaders)
    cHeaders{ii} = sprintf('\\multicolumn{%d}{c}{%s} &',ncols,cHeaders{ii});
end
cHeaders{end} = [cHeaders{end}(1:end-1) '\\'];
cHeaders = [cHeaders; cell(ncols-1,ncols)];
cHeaders = cHeaders(:)';

% Add row headers
rHeaders = {'Excess','\alpha','MKT','SMB','HML'};
rHeaders = [rHeaders; cell(size(rHeaders))];
rHeaders = [{[]}; rHeaders(:)];
rHeaders = strcat(rHeaders, ' &');

out = [rHeaders [cHeaders; out]];

end