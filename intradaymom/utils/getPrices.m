function [st_signal, en_signal, st_hpr, end_hpr] = getPrices(type, s)

nseries                                 = numel(s.permnos);
[st_signal, en_signal, st_hpr, end_hpr] = deal(NaN(1,nseries));
switch type
    case 'taq_exact'

        % Get sampled data
        tmp            = getTaqData([],[],[],[],[],s.datapath,s.mst,false);
        % Price at end of signal
        idx            = serial2hhmmss(tmp.Datetime) == s.END_TIME_SIGNAL;
        [~,col]        = ismember(tmp.Permno(idx), s.permnos);
        en_signal(col) = tmp.Price(idx); 
        % Price at beginning of holding period
        idx            = serial2hhmmss(tmp.Datetime) == s.START_TIME_HPR;
        [~,col]        = ismember(tmp.Permno(idx), s.permnos);
        st_hpr(col)    = tmp.Price(idx);
        
        % Price at beginnig of signal
        [idx,col]      = ismember(s.price_fl.Permno, s.permnos);
        st_signal(col) = s.price_fl.FirstPrice(idx);
        % Price at end of holdong period
        end_hpr(col)   = s.price_fl.LastPrice(idx);
    
    case 'taq_vwap'
    case 'crsp_exact'
    case ''
end
end