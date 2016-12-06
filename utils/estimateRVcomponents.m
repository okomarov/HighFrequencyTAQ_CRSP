function rv = estimateRVcomponents(freq, useon)
if nargin < 1 || isempty(freq),     freq     = 5;    end
if nargin < 2 || isempty(useon),    useon    = true; end

WRITETO   = '.\results\';
FUNNAME   = 'rvcomponents';
PATH2DATA = 'data\TAQ\sampled\5min\nobad_vw';

try
    fname = matname(FUNNAME,freq, useon);
    rv    = loadresults(fname);
catch
    load(fullfile(PATH2DATA,'master'),'-mat','mst');

    % Add overnight return or discard
    if useon
        reton                 = loadresults('return_intraday_overnight','hfandlow\results');
        [idx,pos]             = ismembIdDate(reton.Permno,reton.Date, mst.Permno, mst.Date);
        mst.RetCO(pos(idx),1) = reton.RetCO(idx);
        fun                   = @(x) {mst(x,{'Permno','Date','RetCO'})};
        mst                   = accumarray(mst.File,(1:size(mst))',[], fun);
    end

    % Calculate beta components: sum(r*benchr) and sum(benchr^2)
    fprintf('%s: creating RV components at %d min.\n', mfilename, freq)
    [rv,filename] = Analyze(FUNNAME, [], mst, PATH2DATA,1,8,struct('USE_OVERNIGHT',useon));

    % Rename to append the sampling frequency
    fname       = regexp(filename,'\w+?(?=\.mat)','match','once');
    fname       = [matname(fname,freq, useon),'.mat'];
    newfullname = fullfile(WRITETO, fname);
    movefile(fullfile(WRITETO,filename), newfullname);
end
end

function name = matname(name, freq, useon)
if useon
    useon = 'on';
else
    useon = '';
end
name = sprintf('%s%dm%s%s', name, freq, useon);
end
