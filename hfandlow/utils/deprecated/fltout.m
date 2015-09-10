function [out,mu,mad,lambda] = fltout(dates,price,k,mult)

% FLTOUT filter outliers from data
%
%   FLTOUT(DATA,K,D) Filters outliers contained in DATA (1st col: serial
%                    dates; 2nd col: prices). K are the number of neighbor 
%                    observation around the point examined and D is the 
%                    number of neighbors top/bottom-stripped from the sorted K values
%                    
%   OUT = FILTOUT(...)
%       OUT is a cell array containing the result of H: data -> out
%
% See also:

% Author: Oleg Komarov (oleg.komarov@hotmail.it)
% Tested on R2011a.
% 01 aug 2011 - Created

% NINPUTS
% error(nargchk(4,4,nargin));

% Half neighbours
k  = uint32(k);
k2 = k/2;

% Position of first observations for each day
last = reshape(uint32([0; 
                       find(diff(fix(dates)));
                       size(dates,1)]...
                      ),1,1,[]);
                  
% Begin of the day
Ini = repmat((1:k+1).',1,k2);
Ini(1:k+2:k*k2) = [];
Ini = reshape(Ini, k, k2);
Ini = bsxfun(@plus,Ini,last(1:end-1));

% End of day: rot90(Ini,2)
Fin = bsxfun(@minus, last(2:end), Ini(k:-1:1,k2:-1:1,1))+1;

% Calculate means and std
mu     = zeros(last(end),2);
mad    = zeros(last(end),1);
lambda = zeros(last(end),1);
subin  = [tril(true(k2)); true(k2)]; 
subfi  = [false(k2); tril(true(k2),-1)];

pos1   = [1:k2, k2+2:k+1].';
pos2   = 1:k2;
pos3   = pos2 + k2;

for n = uint16(1:numel(last)-1)
    % Initial and final parts - left/right median
    for r = 1:k2
        mu([last(n)+r last(n+1)-k2+r],:) =...
             [my_median(price(Ini(~subin(:,r),r,n))),...
              my_median(price(Ini( subin(:,r),r,n)));
              my_median(price(Fin(~subfi(:,r),r,n))),...
              my_median(price(Fin( subfi(:,r),r,n)))];
    end
    % Initial and final parts - mean absolute deviation (mad)
    tmp = price(Ini(:,:,n));
    mad(last(n)+1:last(n)+k2) = mean(abs(bsxfun(@minus,tmp,my_median(tmp))));
    tmp = price(Fin(:,:,n));
    mad(last(n+1)-k2+1:last(n+1)) = mean(abs(bsxfun(@minus,tmp,my_median(tmp))));
        
    % Mid part - left/right medians and mad
    neigh      = price(bsxfun(@plus, pos1, uint32(last(n):last(n+1)-k-1)));
    idx        = last(n)+k2+1:last(n+1)-k2;
    mu (idx,:) = [my_median(neigh(pos2,:)); my_median(neigh(pos3,:))].';
    mad(idx)   = mean(abs(bsxfun(@minus,neigh,my_median(neigh))));
    
    % Lambda
    idx = last(n)+1:last(n+1);
    lambda(idx) = mean(abs(diff(price(idx))))*mult;
    
end
mu(isnan(mu)) = 0;
% The filter rule |x - d_mu(left|right)| > max(4*mad,lambda)
lhs = abs(bsxfun(@minus,price,mu));
rhs = max(4*mad,lambda);
out = all(bsxfun(@gt,lhs,rhs),2);
end
function out = my_median(in)
if isempty(in), out = NaN; else out = fast_median(in); end
end