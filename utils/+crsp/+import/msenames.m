function data = msenames(path2zip, outname)
% MSENAMES Import msenames dataset 
%
%     Variable  Type     Label
%     --------  ----     -----    
%     Permno    uint32   PERMNO
%     Namedt    uint32   Names Date (yyyymmdd)
%     Nameendt  uint32   Names Ending Date (yyyymmdd)
%     Shrcd     uint8    Share Code
%     Exchcd    uint8    Exchange Code
%     Siccd     uint16   Standard Industrial Classification Code
%     Ncusip    cellstr  CUSIP historical
%     Ticker    cellstr  Ticker Symbol
%     Comnam    cellstr  Company Name
%     Shrcls    cellstr  Share Class
%     Tsymbol   cellstr  Trading Symbol
%     Naics     uint32   North American Industry Classification System
%     Primexch  cellstr  Primary Exchange
%     Trdstat   cellstr  Trading Status
%     Secstat   cellstr  Security Status
%     Permco    uint32   PERMCO
%     Compno    uint32   Nasdaq Company Number
%     Issuno    uint32   Nasdaq Issue Number
%     Hexcd     uint8    Exchange Code Header
%     Hsiccd    uint16   Standard Industrial Classification Code
%     Cusip     cellstr  CUSIP Header
%
if nargin < 2
    [folder,name] = fileparts(mfilename('fullpath'));
    outname       = fullfile(folder,'..\data',name);
end
fmt  = ['%u32%u32%u32%u8%u8%u16%s%s%s%s '...
        '%s%u32%s%s%s%u32%u32%u32%u8%u16%s'];
labels = {'PERMNO';'Names Date (yyyymmdd)';'Names Ending Date  (yyyymmdd)';'Share Code';'Exchange Code';'Standard Industrial Classification Code';'CUSIP historical';'Ticker Symbol';'Company Name';'Share Class';'Trading Symbol';'North American Industry Classification System';'Primary Exchange';'Trading Status';'Security Status';'PERMCO';'Nasdaq Company Number';'Nasdaq Issue Number';'Exchange Code Header';'Standard Industrial Classification Code';'CUSIP Header'};
data = crsp.import.generic(path2zip, outname, fmt, struct('VariableDescriptions',{labels}));
end
