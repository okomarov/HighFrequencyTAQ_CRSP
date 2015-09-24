function idx = selecttrades(data)

% SELECTTRADES Select invalid TAQ trades
%
%   SELCTTRADES(DATA) where DATA should be a dataset with the 
%                     following varnames (case-sensitive):
%                       - G127_Correction
%                       - Condition
%                       - Price
%
%   NOTE: for the description of the varnames and all possible
%         refer to the .PDFs "Monthly TAQ User Guide v1.3" and 
%         "Daily TAQ Client Specification v1.8"
%
%   IDX = SELECTTRADES(...) Returns the index IDX with the trades
%                           to exclude.
%
% See also:

% Exclude:
% - condition all but 0 (not a “G”, Rule 127, or stopped stock trade) 
%   and 40 (Display Book trade)
% - corrections, i.e. all but 0
% - anomalous sale conditions, i.e. all but 'E' (69) Automatic Execution, 'F' (70) Intermarket Sweep Order, ' ' (32) No Sale Condition required or '@' (64) Regular Trade
% - 0 prices
% - size > 0
% idx = (data.G127_Correction(:,1) ~= 0 & data.G127_Correction(:,1) ~= 40) | data.G127_Correction(:,2) ~= 0 | ~all(data.Condition == ' ' | data.Condition == 'E' | data.Condition == 'F' | data.Condition == '@', 2) | data.Price == 0 | data.Volume <= 0;
idx = all(bsxfun(@ne, data.G127_Correction(:,1), [0,40]),2) | data.G127_Correction(:,2) ~= 0 | (all(bsxfun(@ne, data.Condition(:,1), ' EF@'),2) | all(bsxfun(@ne, data.Condition(:,2), ' EF@'),2)) | data.Price == 0 | data.Volume <= 0;
end