function [data,mst] = queryTAQ(mstentry)
    % Needs fields From, To, File

    % Load data TAQ file
    s          = load(fullfile(cd, '.\data\TAQ', sprintf('T%04d.mat',mstentry.File)));
    % Select data
    data       = s.data(mstentry.From:mstentry.To,:);
    % Select master
    mst        = join(s.mst,mstentry,'Type','Inner','Keys',{'Date','From','To'},'RightVars',[]);
    mst.Ticker = s.ids(mst.Id);
    mst        = mst(:,[end,1:end-1]);
    if nargout < 2
        disp(mst)
    end
end