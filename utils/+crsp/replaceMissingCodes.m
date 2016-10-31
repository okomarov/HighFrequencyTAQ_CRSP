function [out, codes] = replaceMissingCodes(c, val)
% NANMISSINGCODES Sets to NaN missing codes
%
% Missing codes for holding period (RET):
%
%   CRSP code | SAS code | Description
%      -44         .E      No valid comparison for an excess return
%      -55         .D      No listing information
%      -66         .C      No valid previous price or > 10 periods before
%      -77         .B      Not trading on the current exchange
%      -88         .A      No data available to calculate returns
%      -99         .       No valid price (usually suspension or trading on unknown exchange)
%
% Missing codes for delisting return (DLRET):
%
%   CRSP code | SAS code | Description
%      -55        .S       CRSP has no source to establish a value after delisting
%      -66        .T       More than 10 trading periods between a security's last price and its first available price on a new exchange
%      -88        .A       Security is still active
%      -99        .P       Security trades on a new exchange after delisting, but CRSP currently has no sources to gather price information
%
% See for details: https://wrds-web.wharton.upenn.edu/wrds/ds/crsp/stock_a/msf.cfm?navId=128
if nargin < 2
    val = NaN;
end

% Find empty returns
nonempty = ~cellfun('isempty',c);
sz       = size(c);

% Find missing codes in returns and map them out
imisscode           = false(sz);
imisscode(nonempty) = ~cellfun('isempty',regexp(c(nonempty), '^[ABCDEPST\.]','once'));

% Convert into numbers
out      = repmat(val,sz);
idx      = nonempty & ~imisscode;
tmp      = textscan(char(c(idx))', '%9f64');
out(idx) = tmp{1};

if nargout == 2
    codes                = zeros(sz,'uint8');
    [~,codes(imisscode)] = ismember(c(imisscode), {'A';'P';'S';'T'});
end
end