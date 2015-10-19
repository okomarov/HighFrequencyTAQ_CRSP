function [st_signal, en_signal, st_hpr, end_hpr] = getPrices(type, permnos, s)

nseries = numel(permnos);

[st_signal, en_signal, st_hpr, end_hpr] = deal(NaN(1,nseries));

switch type
    case 'taq_exact'

        % Get sampled data
        tmp            = getTaqData([],[],[],[],[],s.datapath,s.mst,false);
        % Price at end of signal
        idx            = serial2hhmmss(tmp.Datetime) == s.END_TIME_SIGNAL;
        [~,col]        = ismember(tmp.Permno(idx), permnos);
        en_signal(col) = tmp.Price(idx); 
        % Price at beginning of holding period
        idx            = serial2hhmmss(tmp.Datetime) == s.START_TIME_HPR;
        [~,col]        = ismember(tmp.Permno(idx), permnos);
        st_hpr(col)    = tmp.Price(idx);
        
        % Price at beginnig of signal
        [idx,col]      = ismember(s.price_fl.Permno, permnos);
        st_signal(col) = s.price_fl.FirstPrice(idx);
        % Price at end of holding period
        end_hpr(col)   = s.price_fl.LastPrice(idx);
    
    case 'taq_vwap'
        [idx,col]      = ismember(s.vwap.Permno, permnos);
        st_signal(col) = s.vwap.T93000(idx);
        en_signal(col) = s.vwap.T120000(idx);
        st_hpr(col)    = s.vwap.T123000(idx);
        end_hpr(col)   = s.vwap.T153000(idx);
    
    case 'taq_exact/vwap'
        % Signal
        [idx,col]      = ismember(s.price_fl.Permno, permnos);
        st_signal(col) = s.price_fl.FirstPrice(idx);
        
        tmp            = getTaqData([],[],[],[],[],s.datapath,s.mst,false);
        idx            = serial2hhmmss(tmp.Datetime) == s.END_TIME_SIGNAL;
        [~,col]        = ismember(tmp.Permno(idx), permnos);
        en_signal(col) = tmp.Price(idx); 

        % HPR
        [idx,col]      = ismember(s.vwap.Permno, permnos);
        st_hpr(col)    = s.vwap.T123000(idx);
        end_hpr(col)   = s.vwap.T153000(idx);
        
    case 'crsp_exact'
        
    
end
end